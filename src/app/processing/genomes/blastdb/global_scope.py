# pylint: disable=
"""
Global genome scope resolution utilities for BLAST database construction.

This module defines the domain-level logic for determining which genome
assemblies should be included in a global BLAST database.
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import Dict, List, Optional


# Public API
def resolve_global_genome_scope(
    *,
    dataset_root: Path,
    genus_list: List[str],
) -> Dict:
    """
    Resolve the global genome inclusion scope for BLAST DB construction.

    This function inspects pre-built genus-level genome datasets and determines
    which genome FASTA files are eligible for inclusion in a global BLAST
    database.

    Inclusion is based on the presence of genome FASTA files under:
        <dataset_root>/<genus>/genomes/fasta/*.fna

    Parameters
    ----------
    dataset_root : pathlib.Path
        Root directory containing registered genus-level genome datasets.
    genus_list : list of str
        List of genus names defining the intended global scope.

    Returns
    -------
    dict
        Structured scope definition with the following schema:
        {
          "included": [
            {
              "genus": <str>,
              "genomes": [
                {
                  "accession": <str>,
                  "fasta_path": <str>
                },
                ...
              ],
              "count": <int>
            },
            ...
          ],
          "excluded": [
            {
              "genus": <str>,
              "reason": <str>
            },
            ...
          ],
          "summary": {
            "total_genera": <int>,
            "included_genera": <int>,
            "excluded_genera": <int>,
            "total_genomes": <int>
          }
        }
    """
    included: List[Dict] = []
    excluded: List[Dict] = []

    total_genomes = 0

    for genus_name in genus_list:
        genus_root = dataset_root / genus_name
        genomes_dir = genus_root / "genomes" / "fasta"

        # Case 1: genus directory does not exist
        if not genus_root.is_dir():
            excluded.append(
                _build_exclusion_record(
                    genus_name,
                    reason="genus_dataset_not_found",
                )
            )
            continue

        # Case 2: genomes/fasta directory missing
        if not genomes_dir.is_dir():
            excluded.append(
                _build_exclusion_record(
                    genus_name,
                    reason="genomes_fasta_directory_not_found",
                )
            )
            continue

        fasta_files = _list_fasta_files(genomes_dir)

        # Case 3: no FASTA files found
        if not fasta_files:
            excluded.append(
                _build_exclusion_record(
                    genus_name,
                    reason="no_genome_fasta_found",
                )
            )
            continue

        genome_records = []
        for fasta_path in fasta_files:
            accession = _parse_accession_from_fasta(fasta_path)
            genome_records.append(
                {
                    "accession": accession,
                    "fasta_path": str(fasta_path),
                }
            )

        count = len(genome_records)
        total_genomes += count

        included.append(
            {
                "genus": genus_name,
                "genomes": genome_records,
                "count": count,
            }
        )

    summary = {
        "total_genera": len(genus_list),
        "included_genera": len(included),
        "excluded_genera": len(excluded),
        "total_genomes": total_genomes,
    }

    return {
        "included": included,
        "excluded": excluded,
        "summary": summary,
    }


# Internal Helpers
def _list_fasta_files(
    genomes_dir: Path,
) -> List[Path]:
    """
    List genome FASTA files under a genomes/fasta directory.

    Parameters
    ----------
    genomes_dir : pathlib.Path
        Path to the genomes/fasta directory.

    Returns
    -------
    list of pathlib.Path
        Sorted list of FASTA file paths.
    """
    return sorted(genomes_dir.glob("*.fna"))


def _parse_accession_from_fasta(
    fasta_path: Path,
) -> Optional[str]:
    """
    Parse genome accession from a FASTA filename.

    This function assumes filenames follow a RefSeq-style convention,
    such as:
        GCF_000016725.1_ASM1672v1_genomic.fna

    Parameters
    ----------
    fasta_path : pathlib.Path
        Path to the FASTA file.

    Returns
    -------
    str or None
        Parsed accession identifier (e.g. "GCF_000016725.1"),
        or None if parsing fails.
    """
    name = fasta_path.stem  # strip file extension, e.g. ".fna"
    parts = name.split("_")

    # Expected at least: <DB>_<ACCESSION>_<...>
    if len(parts) < 2:
        return None

    # e.g. "GCF_000016725.1"
    return f"{parts[0]}_{parts[1]}"


def _build_exclusion_record(
    genus_name: str,
    *,
    reason: str,
) -> Dict:
    """
    Construct an exclusion record for a genus.

    Parameters
    ----------
    genus_name : str
        Scientific genus name.
    reason : str
        Machine-readable exclusion reason.

    Returns
    -------
    dict
        Exclusion record dictionary.
    """
    return {
        "genus": genus_name,
        "reason": reason,
    }
