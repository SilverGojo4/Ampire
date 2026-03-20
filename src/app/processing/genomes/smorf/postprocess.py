"""
Post-processing utilities for smORFinder outputs.

This module provides utilities for cleaning up intermediate files,
normalizing tool-generated output filenames, and exporting smORF
prediction tables in a standardized format used by Ampire datasets.
"""

from __future__ import annotations

# Standard Library Imports
import csv
import glob
import os
import shutil
from pathlib import Path
from typing import Dict, List


# Public API
def cleanup_smorf_tmp(
    *,
    genome_output_dir: Path,
) -> None:
    """
    Remove unnecessary intermediate files from smORFinder tmp directory.

    smORFinder produces several intermediate artifacts that are not
    required for downstream analysis. This function removes selected
    large or redundant files while keeping necessary prediction files.

    Parameters
    ----------
    genome_output_dir : pathlib.Path
        Path to the genome-level smORFinder output directory.

    Returns
    -------
    None
    """

    tmp_dir = genome_output_dir / "tmp"

    if not tmp_dir.exists():
        return

    remove_targets = [
        "prodigal.faa",
        "prodigal.ffn",
        "prodigal.gff",
    ]

    for filename in remove_targets:

        file_path = tmp_dir / filename

        if file_path.exists():
            file_path.unlink()


def rename_smorf_outputs(
    *,
    genome_output_dir: Path,
) -> None:
    """
    Normalize smORFinder output filenames.

    smORFinder generates filenames containing unpredictable prefixes
    derived from genome identifiers. This function renames these files
    into stable, canonical names used by Ampire pipelines.

    Parameters
    ----------
    genome_output_dir : pathlib.Path
        Path to the genome-level smORFinder output directory.

    Returns
    -------
    None
    """
    # Define renaming rules
    rename_map = {
        ".faa": "peptides.faa",
        ".ffn": "nucleotides.ffn",
        ".gff": "annotations.gff",
        ".tsv": "predictions.tsv",
    }

    # Apply renaming
    for ext, new_name in rename_map.items():
        matches = glob.glob(os.path.join(genome_output_dir, f"*{ext}"))
        if not matches:
            continue
        old_path = matches[0]
        new_path = os.path.join(genome_output_dir, new_name)

        # Avoid overwriting if target already exists
        if os.path.exists(new_path):
            continue

        shutil.move(old_path, new_path)


def save_smorf_table(
    *,
    records: List[Dict[str, str]],
    output_path: Path,
) -> None:
    """
    Save smORF prediction records as a TSV table.

    Parameters
    ----------
    records : list[dict[str, str]]
        List of smORF prediction records.
    output_path : pathlib.Path
        Output TSV file path.

    Returns
    -------
    None
    """

    if not records:
        return

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write header and rows
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def convert_predictions_tsv_to_csv(
    *,
    genome_output_dir: Path,
) -> None:
    """
    Convert predictions.tsv into predictions.csv.

    Some downstream workflows prefer CSV format. This function
    converts the TSV prediction table generated during assembly
    into a CSV file while preserving the original content.

    Parameters
    ----------
    genome_output_dir : pathlib.Path
        Path to the genome-level smORF output directory.

    Returns
    -------
    None
    """

    tsv_path = genome_output_dir / "predictions.tsv"
    csv_path = genome_output_dir / "predictions.csv"

    if os.path.exists(tsv_path):
        # Read TSV manually
        with open(tsv_path, "r", newline="", encoding="utf-8") as f_in:
            reader = csv.reader(f_in, delimiter="\t")
            rows = list(reader)

        # Write as CSV (comma-separated)
        with open(csv_path, "w", newline="", encoding="utf-8") as f_out:
            writer = csv.writer(f_out, delimiter=",")
            writer.writerows(rows)

        # Optionally remove old .tsv to avoid redundancy
        os.remove(tsv_path)
