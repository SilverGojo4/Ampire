# pylint: disable=too-many-arguments, too-many-locals
"""
Genome download utilities for genus-level datasets.

This module implements the production-stage logic for downloading all
available genome assemblies for a given genus using NCBI resources.
"""

from __future__ import annotations

# Standard Library Imports
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import List

# Third-Party Library Imports
import pandas as pd

# App Imports
from app.utils.fs_utils import decompress_gzip_files, remove_directory_tree


# Public API
def download_genomes_from_genus(
    *,
    genus_name: str,
    species_table: pd.DataFrame,
    data_root: Path,
    assembly_level: str = "complete,chromosome",
    refseq_categories: str = "all",
    source: str = "refseq",
    parallel: int = 10,
    batch_size: int = 3000,
    overwrite: bool = False,
) -> tuple[pd.DataFrame, Path]:
    """
    Download genome assemblies for a given genus.

    This function retrieves all available genome assemblies
    corresponding to the resolved species of a genus using
    `ncbi-genome-download`, and stores the downloaded FASTA
    files within the Ampire data hierarchy.

    Parameters
    ----------
    genus_name : str
        Scientific genus name (e.g., "Escherichia").
    species_table : pd.DataFrame
        Species-level taxonomy table returned by
        `resolve_species_from_genus`.
    data_root : pathlib.Path
        Root directory for data storage (e.g., data/interim).
    assembly_level : str, optional
        Assembly level filters passed to
        `ncbi-genome-download` (e.g., "complete,chromosome").
    refseq_categories : str, optional
        RefSeq category filters applied during genome retrieval.
    source : str, optional
        Genome source database (e.g., "refseq" or "genbank").
    parallel : int, optional
        Number of parallel download workers.
    batch_size : int, optional
        Maximum number of assemblies processed per download batch.
    overwrite : bool, optional
        Whether to remove existing genome FASTA files before
        downloading. If False, existing files are preserved.

    Returns
    -------
    Tuple[pd.DataFrame, pathlib.Path]
        A tuple of (assemblies_df, genomes_root).
        - If genome assemblies are found, assemblies_df contains
        the normalized metadata table.
        - If no genome assemblies match the retrieval criteria,
        assemblies_df is empty.

    Raises
    ------
    ValueError
        If the input species table is empty or contains no valid
        species-level taxonomic identifiers. This indicates an
        invalid or incomplete taxonomy resolution prior to genome
        retrieval.
    RuntimeError
        If genome download fails due to unexpected errors during
        execution of `ncbi-genome-download`, such as command
        execution failures, I/O errors, or metadata parsing issues.
        These errors are propagated from the underlying download
        routine without being caught or suppressed.
    """

    # Extract species taxids
    if species_table.empty:
        raise ValueError("No species available for genome download.")

    species_taxids: List[int] = (
        species_table["species_taxid"].dropna().astype(int).tolist()
    )
    if not species_taxids:
        raise ValueError("No valid species_taxid found in species_table.")

    # Construct output directories
    genus_root = data_root / genus_name
    genomes_root = genus_root / "genomes"
    fasta_dir = genomes_root / "fasta"

    # Decide whether to execute download
    if overwrite and fasta_dir.exists():
        remove_directory_tree(str(fasta_dir), include_self=True)

    # Ensure fresh fasta directory
    fasta_dir.mkdir(parents=True, exist_ok=True)

    # Execute download (real download)
    raw_metadata = _run_ncbi_genome_download(
        species_taxids=species_taxids,
        output_fasta_dir=fasta_dir,
        assembly_level=assembly_level,
        refseq_categories=refseq_categories,
        source=source,
        parallel=parallel,
        batch_size=batch_size,
    )

    if raw_metadata.empty:
        return pd.DataFrame(), genomes_root

    # Normalize + persist canonical assemblies.csv
    assemblies_df = _normalize_assemblies_metadata(
        df=raw_metadata,
        genus_name=genus_name,
        filters_applied=_format_filters(
            assembly_level=assembly_level,
            refseq_categories=refseq_categories,
            source=source,
        ),
    )

    return assemblies_df, genomes_root


# Internal Helpers
def _run_ncbi_genome_download(
    *,
    species_taxids: List[int],
    output_fasta_dir: Path,
    assembly_level: str,
    refseq_categories: str,
    source: str,
    parallel: int,
    batch_size: int,
) -> pd.DataFrame:
    """
    Execute genome downloads using ncbi-genome-download.

    Parameters
    ----------
    species_taxids : list[int]
        List of NCBI species-level taxonomic identifiers used
        as input for genome retrieval.
    output_fasta_dir : pathlib.Path
        Target directory where genome FASTA files and temporary
        download artifacts are written.
    assembly_level : str
        Assembly level filters passed to ncbi-genome-download
        (e.g., "complete,chromosome").
    refseq_categories : str
        RefSeq category filters applied during genome retrieval.
    source : str
        Genome source database (e.g., "refseq" or "genbank").
    parallel : int
        Number of parallel download workers.
    batch_size : int
        Number of species taxids processed per download batch.

    Returns
    -------
    pd.DataFrame
        Merged genome assembly metadata table produced by
        ncbi-genome-download. An empty DataFrame is returned
        if no assemblies are retrieved.

    Raises
    ------
    RuntimeError
        If ncbi-genome-download terminates with an unexpected
        error during execution.
    """
    if not species_taxids:
        return pd.DataFrame()

    assembly_level = assembly_level.replace(" ", "")
    refseq_categories = refseq_categories.replace(" ", "")

    required_cols = [
        "assembly_accession",
        "refseq_category",
        "relation_to_type_material",
        "taxid",
        "species_taxid",
        "organism_name",
        "infraspecific_name",
        "assembly_level",
        "ftp_path",
        "local_filename",
    ]

    batch_metadata_files: List[Path] = []

    total = len(species_taxids)
    tmp_dir = output_fasta_dir / ".tmp"
    tmp_dir.mkdir(exist_ok=True)
    for batch_index in range(0, total, batch_size):
        batch = species_taxids[batch_index : batch_index + batch_size]
        taxid_str = ",".join(map(str, batch))

        batch_no = batch_index // batch_size + 1
        batch_metadata_path = tmp_dir / f"genomes_batch_{batch_no}.tsv"

        cmd: List[str] = [
            "ncbi-genome-download",
            "bacteria",
            "--species-taxids",
            taxid_str,
            "--refseq-categories",
            refseq_categories,
            "--assembly-level",
            assembly_level,
            "--formats",
            "fasta",
            "--output-folder",
            str(output_fasta_dir),
            "--metadata-table",
            str(batch_metadata_path),
            "--flat-output",
            "--parallel",
            str(parallel),
        ]

        _add_ncbi_source_flag(cmd, source)

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        stderr = (result.stderr or "").strip()
        if result.returncode != 0 and "No downloads matched your filter" not in stderr:
            raise RuntimeError(
                "ncbi-genome-download failed.\n"
                f"Command: {' '.join(cmd)}\n"
                f"STDERR: {stderr}"
            )

        if batch_metadata_path.exists() and batch_metadata_path.stat().st_size > 0:
            batch_metadata_files.append(batch_metadata_path)

    if not batch_metadata_files:
        return pd.DataFrame()

    # Merge metadata
    dfs = [_read_ncbi_metadata_table(p, required_cols) for p in batch_metadata_files]

    df_meta = (
        pd.concat(dfs, ignore_index=True)
        .drop_duplicates(subset=["assembly_accession"])
        .reset_index(drop=True)
    )

    # Cleanup batch metadata TSVs
    for p in batch_metadata_files:
        p.unlink(missing_ok=True)

    # Clean local_filename
    df_meta["local_filename"] = (
        df_meta["local_filename"]
        .astype(str)
        .str.replace(r"^\./", "", regex=True)
        .str.replace(r"\.gz$", "", regex=True)
    )

    # Remove human_readable directory if present
    hr_dir = output_fasta_dir / "human_readable"
    if hr_dir.exists():
        remove_directory_tree(str(hr_dir), include_self=True)

    # Decompress all .gz genome files
    decompress_gzip_files(str(output_fasta_dir), recursive=True)
    if tmp_dir.exists():
        remove_directory_tree(str(tmp_dir), include_self=True)

    return df_meta


def _normalize_assemblies_metadata(
    *,
    df: pd.DataFrame,
    genus_name: str,
    filters_applied: str,
) -> pd.DataFrame:
    """
    Normalize raw genome metadata into the assemblies dataset schema.

    Parameters
    ----------
    df : pd.DataFrame
        Raw genome assembly metadata table returned by
        genome download utilities (e.g., ncbi-genome-download).
    genus_name : str
        Scientific genus name associated with the assemblies.
    filters_applied : str
        String representation of filters applied during genome
        retrieval, recorded for dataset provenance.

    Returns
    -------
    pd.DataFrame
        Normalized assembly metadata table conforming to the
        assemblies.csv schema.
    """
    df = df.copy()

    df["genus_name"] = genus_name
    df["filters_applied"] = filters_applied
    df["retrieved_at"] = datetime.now(timezone.utc).isoformat()

    columns = [
        "assembly_accession",
        "taxid",
        "species_taxid",
        "organism_name",
        "infraspecific_name",
        "genus_name",
        "assembly_level",
        "refseq_category",
        "relation_to_type_material",
        "source",
        "ftp_path",
        "local_filename",
        "retrieved_at",
        "filters_applied",
    ]

    return df.reindex(columns=columns).reset_index(drop=True)


def _read_ncbi_metadata_table(
    file_path: Path,
    required_columns: List[str],
) -> pd.DataFrame:
    """
    Read and validate a metadata table produced by ncbi-genome-download.

    Parameters
    ----------
    file_path : pathlib.Path
        Path to the metadata TSV file generated by
        ncbi-genome-download.
    required_columns : list[str]
        List of column names that must be present in the metadata
        table for successful parsing.

    Returns
    -------
    pd.DataFrame
        Metadata table restricted to the required columns.
        A copy of the data is returned to prevent unintended
        side effects.

    Raises
    ------
    ValueError
        If the metadata table does not contain all required
        schema columns expected by downstream processors.
    """
    df = pd.read_csv(file_path, sep="\t")

    missing = [c for c in required_columns if c not in df.columns]
    if missing:
        raise ValueError(
            "Unexpected metadata schema from ncbi-genome-download. "
            f"Missing columns: {missing}"
        )

    return df[required_columns].copy()


def _format_filters(
    *,
    assembly_level: str,
    refseq_categories: str,
    source: str,
) -> str:
    """
    Format genome retrieval filters as a provenance string.

    Parameters
    ----------
    assembly_level : str
        Assembly level filters applied during genome retrieval.
    refseq_categories : str
        RefSeq category filters applied during genome retrieval.
    source : str
        Genome source database (e.g., "refseq" or "genbank").

    Returns
    -------
    str
        Semicolon-delimited filter string recorded in
        assemblies metadata and dataset manifests.
    """
    return (
        f"assembly_level={assembly_level};"
        f"refseq_categories={refseq_categories};"
        f"source={source}"
    )


def _add_ncbi_source_flag(
    cmd: List[str],
    source: str,
) -> None:
    """
    Append genome source selection flags to an ncbi-genome-download command.

    Parameters
    ----------
    cmd : list[str]
        Command argument list used to invoke ncbi-genome-download.
        The list is modified in place.
    source : str
        Genome source database to retrieve assemblies from.
        Supported values are "refseq", "genbank", or "all".
        When set to "all", no source flag is appended.

    Returns
    -------
    None
        This function modifies the command argument list in place.

    Raises
    ------
    ValueError
        If an unsupported genome source is specified.
    """
    if source == "all":
        return

    if source not in {"refseq", "genbank"}:
        raise ValueError(f"Unsupported source: '{source}'")

    # Use short flag for compatibility
    cmd.extend(["-s", source])
