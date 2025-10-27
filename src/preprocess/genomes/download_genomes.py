# pylint: disable=import-error, wrong-import-position, too-many-arguments, too-many-positional-arguments, too-many-locals
"""
Download bacterial genomes from NCBI RefSeq based on genus-level taxonomy.
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import subprocess
import sys
import time
from typing import List, Tuple

# ============================== Third-Party Library Imports ==============================
import pandas as pd
from ete3 import NCBITaxa

# ============================== Project Root Path Setup ==============================
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
LOGGING_PATH = os.path.join(BASE_PATH, "src/utils/logging_toolkit/src/python")

# Fix Python Path
if BASE_PATH not in sys.path:
    sys.path.append(BASE_PATH)
if LOGGING_PATH not in sys.path:
    sys.path.append(LOGGING_PATH)

# ============================== Project-Specific Imports ==============================
# Submodule imports (external tools integrated into project)
from setup_logging import CustomLogger

# Local Ampire utility modules
from src.utils.io_utils import (
    decompress_all_gz_files,
    directory_exists,
    load_dataframe_by_columns,
    remove_directory_recursively,
)
from src.utils.log_utils import (
    get_pipeline_completion_message,
    get_task_completion_message,
)


# ============================== Custom Functions ==============================
def get_species_from_genus(
    genus_name: str,
    output_path: str,
    logger: CustomLogger,
) -> Tuple[pd.DataFrame, bool]:
    """
    Retrieve all species belonging to a specified bacterial genus using ETE3.

    Parameters
    ----------
    genus_name : str
        Name of the genus to query (e.g., "Escherichia").
    output_path : str
        Path to save the resulting DataFrame as a CSV file.
    logger : CustomLogger
        Logger instance for progress tracking.

    Returns
    -------
    (pd.DataFrame, bool)
        A tuple (df_species, success):
            - df_species: DataFrame containing ['Species', 'TaxID'] (empty if not found)
            - success: True if genus found and species retrieved, False otherwise
    """
    logger.info("/ Task: Retrieve species taxonomy from genus")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Ensure output directory exists
        metadata_output_dir = os.path.dirname(output_path)
        os.makedirs(metadata_output_dir, exist_ok=True)

        # Initialize NCBI taxonomy database
        ncbi = NCBITaxa()

        # Retrieve genus TaxID
        name_map = ncbi.get_name_translator([genus_name])
        if not name_map or genus_name not in name_map:
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    f"Genus '{genus_name}' not found in NCBI taxonomy. "
                    f"Skipping this genus."
                ),
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            df_species = pd.DataFrame(columns=["Species", "TaxID"])

        else:
            # Retrieve genus taxid
            genus_id = name_map[genus_name][0]

            # Fetch all descendant taxa under the genus
            descendants = ncbi.get_descendant_taxa(genus_id, intermediate_nodes=True)
            ranks = ncbi.get_rank(descendants)
            names = ncbi.get_taxid_translator(descendants)

            # Keep only "species" level taxa
            records = [
                {"Species": names[t], "TaxID": t}
                for t in descendants
                if ranks.get(t) == "species"
            ]

            # Convert to DataFrame and sort alphabetically
            df_species = (
                pd.DataFrame(records).sort_values("Species").reset_index(drop=True)
            )
            n_species = df_species.shape[0]

            # Log summary
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    f"[ Taxonomy summary ]\n"
                    f"▸ Genus   : '{genus_name}'\n"
                    f"▸ Species : {n_species}"
                ),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

            # Save to CSV
            df_species.to_csv(output_path, index=False)
            logger.log_with_borders(
                level=logging.INFO,
                message=f"Saved:\n'{genus_name}/metadata/species.csv'",
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Log completion
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Determine success flag
        success = not df_species.empty
        return df_species, success

    except Exception:
        logger.exception("Unexpected error in 'get_species_from_genus()'")
        raise


def ncbi_genome_download(
    species_taxids: List[int],
    output_dir: str,
    logger: CustomLogger,
    assembly_level: str = "complete,chromosome",
    refseq_categories: str = "all",
    formats: str = "fasta",
    parallel: int = 10,
    verbose: bool = True,
    batch_size: int = 3000,
) -> None:
    """
    Run NCBI Genome Download (bacteria) for a given list of species TaxIDs.

    Parameters
    ----------
    species_taxids : List[int]
        List of NCBI TaxIDs corresponding to species to download.
    output_dir : str
        Directory where genome files and metadata will be stored.
    logger : CustomLogger
        Project-specific logger for structured logging.
    assembly_level : str, optional
        Genome assembly level to include (default: "complete,chromosome").
    refseq_categories : str, optional
        RefSeq categories to include (default: "all").
    formats : str, optional
        File formats to download (default: "fasta").
    parallel : int, optional
        Number of parallel download threads (default: 10).
    verbose : bool, optional
        Whether to enable verbose mode (default: True).
    batch_size : int, optional
        Number of TaxIDs to include per batch download. This prevents command
        invocation from exceeding system argument-size limits. (default: 3000)

    Returns
    -------
    None
    """
    logger.info("/ Task: Run 'ncbi-genome-download' for bacterial genomes")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Ensure output directory exists
        fasta_output_dir = os.path.join(output_dir, "fasta")
        metadata_output_dir = os.path.join(output_dir, "metadata")
        for sub_output_dir in [fasta_output_dir, metadata_output_dir]:
            os.makedirs(sub_output_dir, exist_ok=True)

        # Prepare output metadata path
        metadata_path = os.path.join(metadata_output_dir, "genomes.tsv")

        # Log parameters and command
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ NCBI Genome Download Parameters ]\n"
                f"▸ Assembly level   : {assembly_level}\n"
                f"▸ RefSeq category  : {refseq_categories}\n"
                f"▸ Formats          : {formats}\n"
                f"▸ Parallel threads : {parallel}\n"
                f"▸ Verbose mode     : {verbose}\n"
                f"▸ Batch size       : {batch_size if batch_size else 'Disabled'}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        assembly_level = assembly_level.replace(" ", "")

        # Always run in batch mode — prevents argument too long
        total = len(species_taxids)
        batch_metadata_files = []
        for batch_index in range(0, total, batch_size):

            # Convert TaxID list to comma-separated string
            batch = species_taxids[batch_index : batch_index + batch_size]
            taxid_str = ",".join(map(str, batch))

            batch_metadata_path = os.path.join(
                metadata_output_dir, f"genomes_batch_{batch_index//batch_size+1}.tsv"
            )

            # Construct command arguments
            cmd = [
                "ncbi-genome-download",
                "bacteria",
                "--species-taxids",
                taxid_str,
                "--refseq-categories",
                refseq_categories,
                "--assembly-level",
                assembly_level,
                "--formats",
                formats,
                "--output-folder",
                fasta_output_dir,
                "--metadata-table",
                batch_metadata_path,
                "--human-readable",
                "--flat-output",
                "--parallel",
                str(parallel),
            ]

            if verbose:
                cmd.append("--verbose")

            # Execute the command
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if (
                result.returncode != 0
                and "No downloads matched your filter" not in result.stderr
            ):
                raise RuntimeError(result.stderr)

            # Check metadata existence only when not a fatal failure
            if (
                os.path.exists(batch_metadata_path)
                and os.path.getsize(batch_metadata_path) > 0
            ):
                batch_metadata_files.append(batch_metadata_path)

        if len(batch_metadata_files) == 0:
            # No metadata generated across all batches
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    "No genomes matched the given filters.\n"
                    "This genus may not have entries in RefSeq or GenBank."
                ),
                border="|",
                length=120,
            )

        else:

            if len(batch_metadata_files) == 1:

                # Only one metadata file → rename to genomes.tsv directly
                single_meta_path = batch_metadata_files[0]
                df_meta = load_dataframe_by_columns(file_path=single_meta_path)
                df_meta.to_csv(metadata_path, sep="\t", index=False)
                os.remove(single_meta_path)

            else:

                # Multiple metadata batches → concatenate and deduplicate
                df_list = []
                for meta_file in batch_metadata_files:
                    df_list.append(load_dataframe_by_columns(file_path=meta_file))

                # Concatenate and drop duplicates
                df_meta = pd.concat(df_list, ignore_index=True).drop_duplicates()

                # Output final merged metadata to TSV
                df_meta.to_csv(metadata_path, sep="\t", index=False)

                # Remove intermediate batch files
                for meta_file in batch_metadata_files:
                    os.remove(meta_file)

            # Verify metadata file existence
            if not os.path.exists(metadata_path):
                raise RuntimeError(f"Metadata file not found: '{metadata_path}'")

            # Remove "human_readable" folder if exists
            hr_dir = os.path.join(fasta_output_dir, "human_readable")
            remove_directory_recursively(dir_path=hr_dir)

            # Decompress all downloaded .gz genome files
            decompress_all_gz_files(input_dir=fasta_output_dir, recursive=True)

            # File exists — read TSV
            required_cols = [
                "assembly_accession",
                "refseq_category",
                "relation_to_type_material",
                "taxid",
                "species_taxid",
                "organism_name",
                "infraspecific_name",
                "assembly_level",
                "seq_rel_date",
                "ftp_path",
                "local_filename",
            ]
            df_meta = load_dataframe_by_columns(
                file_path=metadata_path, required_columns=required_cols
            )
            n_records = df_meta.shape[0]

            # Clean 'local_filename' column — remove leading './' and trailing '.gz'
            df_meta["local_filename"] = (
                df_meta["local_filename"]
                .str.replace(r"^\./", "", regex=True)
                .str.replace(r"\.gz$", "", regex=True)
            )

            # Convert TSV to CSV
            csv_path = metadata_path.replace(".tsv", ".csv")
            parts = csv_path.strip("/").split("/")
            rel_path = "/".join(parts[-3:])
            df_meta.to_csv(csv_path, index=False)
            os.remove(metadata_path)

            # Log metadata summary
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    "[ Metadata Summary ]\n" f"▸ Records downloaded : {n_records}\n"
                ),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            logger.log_with_borders(
                level=logging.INFO,
                message=f"Saved:\n'{rel_path}'",
                border="|",
                length=120,
            )

        # Log completion
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error during 'run_ncbi_genome_download()'.")
        raise


# ============================== Pipeline Entry Point ==============================
def run_download_genomes_by_genus(
    base_path: str, logger: CustomLogger, **kwargs
) -> None:
    """
    Download all bacterial genomes from NCBI RefSeq for a specified genus.

    Parameters
    ----------
    base_path : str
        Project root path.
    logger : CustomLogger
        Logger instance for tracking progress.
    **kwargs : dict
        Optional overrides for input/output file paths.

    Returns
    -------
    None
    """
    # Start timing
    start_time = time.time()

    try:
        # -------------------- Retrieve input parameters (CLI or default) --------------------
        genus_name = kwargs.get("genomes_genus_name", "Escherichia")
        output_dir = kwargs.get(
            "genomes_output_dir",
            os.path.join(base_path, "Escherichia"),
        )

        # Ensure output directory exists
        if not directory_exists(dir_path=output_dir):
            os.makedirs(name=output_dir, exist_ok=True)

        # Retrieve species under the genus
        species_csv = os.path.join(output_dir, "metadata/species.csv")
        df_species, success = get_species_from_genus(genus_name, species_csv, logger)
        logger.add_spacer(level=logging.INFO, lines=1)

        # Download genomes
        if success:
            species_taxids = df_species["TaxID"].tolist()
            ncbi_genome_download(species_taxids, output_dir, logger)
            logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_download_genomes_by_genus()'")
        raise

    # Final summary block
    logger.info("[ 'Pipeline Execution Summary' ]")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="=")
    logger.log_with_borders(
        level=logging.INFO,
        message="\n".join(get_pipeline_completion_message(start_time)),
        border="║",
        length=120,
    )
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="=")
