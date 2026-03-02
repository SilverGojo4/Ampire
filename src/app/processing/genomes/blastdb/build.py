# pylint: disable=too-many-locals
"""
Global BLAST database construction utilities.

This module defines the domain-level logic for constructing a global
BLAST database from a resolved genome inclusion scope.
"""

from __future__ import annotations

# Standard Library Imports
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List


# Public API
def build_global_blast_db(
    *,
    scope: Dict,
    output_dir: Path,
    db_name: str = "ampire_genomes",
    molecule_type: str = "nucl",
    overwrite: bool = False,
) -> Dict:
    """
    Build a global BLAST database from a resolved genome scope.

    This function materializes a global genome FASTA with normalized,
    globally unique sequence identifiers and constructs a BLAST database
    using `makeblastdb`.

    Parameters
    ----------
    scope : dict
        Global genome scope definition returned by
        `resolve_global_genome_scope()`.
    output_dir : pathlib.Path
        Directory where BLAST database files will be created.
    db_name : str, optional
        Base name of the BLAST database.
    molecule_type : {"nucl", "prot"}, optional
        Molecule type for BLAST database construction.
    overwrite : bool, optional
        Whether to overwrite existing BLAST database files.

    Returns
    -------
    dict
        Global BLAST database build summary with the following structure:
        {
            "blastdb": {
                "name": str,
                "type": str,
                "scope": "global",
                "status": "built" | "reused",
            },
            "source": {
                "genera_included": int,
                "genomes_included": int,
                "total_bases": int,
                "excluded_genera": list[str],
            },
            "artifacts": {
                "db_base": str,
                "fasta_path": str,
            },
            "provenance": {
                "tool": "makeblastdb",
                "molecule_type": str,
            },
        }

    Raises
    ------
    RuntimeError
        If BLAST database construction fails.
    """

    # Resolve and prepare output paths
    output_dir = output_dir.resolve()
    blastdb_dir = output_dir / db_name

    # Overwrite mode: remove entire blastdb directory to avoid stale fragments
    if overwrite and blastdb_dir.exists():
        shutil.rmtree(blastdb_dir)

    # Ensure blastdb directory exists
    blastdb_dir.mkdir(parents=True, exist_ok=True)
    db_base = blastdb_dir / db_name
    concat_fasta = blastdb_dir / f"{db_name}.fasta"

    # Determine reuse vs rebuild
    exists = _blast_db_exists(db_base=db_base, molecule_type=molecule_type)

    if exists and not overwrite:
        status = "reused"
    else:
        status = "built"

        # 1. Collect FASTA files and genus mapping from scope
        fasta_files: List[Path] = []
        genus_map: Dict[Path, str] = {}

        included = scope.get("included", [])

        for genus_entry in included:
            genus = genus_entry.get("genus")
            for genome in genus_entry.get("genomes", []):
                fasta_path = genome.get("fasta_path")
                if fasta_path:
                    path = Path(fasta_path)
                    fasta_files.append(path)
                    genus_map[path] = genus

        if not fasta_files:
            raise RuntimeError(
                "Global genome scope contains no FASTA files; "
                "BLAST database construction aborted."
            )

        # 2. Write concatenated FASTA with normalized headers
        _write_concatenated_fasta(
            fasta_files=fasta_files,
            output_fasta=concat_fasta,
            genus_map=genus_map,
        )

        # 3. Build BLAST database
        _run_makeblastdb(
            fasta_path=concat_fasta,
            db_base=db_base,
            molecule_type=molecule_type,
        )

    # 4. Compute source summary statistics
    genera_included = len(scope.get("included", []))
    excluded_genera = [g.get("genus") for g in scope.get("excluded", [])]
    genomes_included = 0
    total_bases = 0

    for genus_entry in scope.get("included", []):
        for genome in genus_entry.get("genomes", []):
            genomes_included += 1
            total_bases += genome.get("length", 0)

    # 5. Return structured build summary
    return {
        "blastdb": {
            "name": db_name,
            "type": molecule_type,
            "scope": "global",
            "status": status,
        },
        "source": {
            "genera_included": genera_included,
            "genomes_included": genomes_included,
            "total_bases": total_bases,
            "excluded_genera": excluded_genera,
        },
        "artifacts": {
            "db_base": str(db_base),
            "fasta_path": str(concat_fasta),
        },
        "provenance": {
            "tool": "makeblastdb",
            "molecule_type": molecule_type,
        },
    }


# Internal helpers
def _write_concatenated_fasta(
    *,
    fasta_files: Iterable[Path],
    output_fasta: Path,
    genus_map: Dict[Path, str],
) -> None:
    """
    Concatenate multiple FASTA files into a single FASTA file with
    normalized headers.

    Parameters
    ----------
    fasta_files : iterable of pathlib.Path
        Input genome FASTA file paths.
    output_fasta : pathlib.Path
        Output FASTA file to be written.
    genus_map : dict[pathlib.Path, str]
        Mapping from FASTA file path to genus name.

    Raises
    ------
    RuntimeError
        If a duplicated normalized sequence identifier is encountered
        during FASTA concatenation.
    """
    seen_seqids: set[str] = set()

    with output_fasta.open("w") as out_handle:
        for fasta_path in fasta_files:
            genus = genus_map.get(fasta_path)
            if genus is None:
                raise RuntimeError(
                    f"Genus mapping not found for FASTA file: '{fasta_path}'"
                )

            with fasta_path.open("r") as in_handle:
                current_seqid: str | None = None

                for line in in_handle:
                    if line.startswith(">"):
                        raw_header = line[1:].strip()

                        normalized_seqid, normalized_header = _normalize_fasta_header(
                            raw_header=raw_header,
                            genus=genus,
                        )

                        # Enforce global uniqueness
                        if normalized_seqid in seen_seqids:
                            raise RuntimeError(
                                "Duplicated sequence identifier encountered during "
                                "global FASTA construction: "
                                f"'{normalized_seqid}' (source: '{fasta_path}')"
                            )

                        seen_seqids.add(normalized_seqid)
                        current_seqid = normalized_seqid

                        out_handle.write(f">{normalized_header}\n")

                    else:
                        # Sequence line
                        if current_seqid is None:
                            # Malformed FASTA: sequence before header
                            raise RuntimeError(
                                f"Malformed FASTA file encountered: '{fasta_path}'"
                            )
                        out_handle.write(line)


def _normalize_fasta_header(
    *,
    raw_header: str,
    genus: str,
) -> tuple[str, str]:
    """
    Normalize a FASTA header according to the Ampire global
    BLAST database header contract.

    Header contract (v0.1):
        >{accession}|genus={genus} {optional_original_description}

    Parameters
    ----------
    raw_header : str
        Original FASTA header line without the leading '>'.
    genus : str
        Taxonomic genus name associated with the sequence.

    Returns
    -------
    tuple[str, str]
        (normalized_seqid, normalized_header)
        - normalized_seqid: machine-readable unique identifier
        - normalized_header: full header string (without leading '>')

    Raises
    ------
    RuntimeError
        If the accession cannot be parsed from the raw header.
    """
    raw_header = raw_header.strip()
    if not raw_header:
        raise RuntimeError("Encountered empty FASTA header.")

    # Accession = first whitespace-delimited token
    parts = raw_header.split(maxsplit=1)
    accession = parts[0]

    if not accession:
        raise RuntimeError(
            f"Failed to parse accession from FASTA header: '{raw_header}'"
        )

    normalized_seqid = f"{accession}|genus={genus}"

    # Preserve original description if present
    if len(parts) > 1:
        description = parts[1]
        normalized_header = f"{normalized_seqid} {description}"
    else:
        normalized_header = normalized_seqid

    return normalized_seqid, normalized_header


def _run_makeblastdb(
    *,
    fasta_path: Path,
    db_base: Path,
    molecule_type: str,
) -> None:
    """
    Run makeblastdb to construct a BLAST database.

    Parameters
    ----------
    fasta_path : pathlib.Path
        Path to concatenated FASTA file.
    db_base : pathlib.Path
        Base path for BLAST database output.
    molecule_type : {"nucl", "prot"}
        Molecule type for BLAST database.
    """
    cmd = [
        "makeblastdb",
        "-in",
        str(fasta_path),
        "-dbtype",
        molecule_type,
        "-out",
        str(db_base),
        "-parse_seqids",
    ]

    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "makeblastdb failed:\n" f"STDOUT:\n{exc.stdout}\n" f"STDERR:\n{exc.stderr}"
        ) from exc


def _blast_db_exists(
    *,
    db_base: Path,
    molecule_type: str,
) -> bool:
    """
    Check whether a BLAST database exists (strong check).

    This function validates BLAST database existence by checking
    the required index files generated by `makeblastdb`.

    Required file sets:
    - nucl: .nin, .nsq, .nhr
    - prot: .pin, .psq, .phr

    Parameters
    ----------
    db_base : pathlib.Path
        Base path of BLAST database.
    molecule_type : {"nucl", "prot"}
        Molecule type for BLAST database.

    Returns
    -------
    bool
        True if required BLAST database files exist.

    Raises
    ------
    ValueError
        If `molecule_type` is invalid.
    """
    if molecule_type == "nucl":
        required_suffixes = [".nin", ".nsq", ".nhr"]
        return all(db_base.with_suffix(sfx).is_file() for sfx in required_suffixes)

    if molecule_type == "prot":
        required_suffixes = [".pin", ".psq", ".phr"]
        return all(db_base.with_suffix(sfx).is_file() for sfx in required_suffixes)

    raise ValueError(f"Invalid molecule_type: '{molecule_type}'")
