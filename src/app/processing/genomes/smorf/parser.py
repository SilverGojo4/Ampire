# pylint: disable=
"""
Parsing utilities for smORFinder outputs.

This module provides robust parsers for smORFinder intermediate
files and assembles the final smORF prediction records used by
Ampire pipelines.

The parsing logic intentionally avoids strict schema assumptions
to maintain compatibility across smORFinder versions.
"""

from __future__ import annotations

# Standard Library Imports
import csv
from pathlib import Path
from typing import Dict, List


# Public API
def read_tsv(
    path: Path,
    delimiter: str = "\t",
    comment: str = "#",
) -> List[Dict[str, str]]:
    """
    Read a delimited text file into a list of dictionaries.

    This function intentionally avoids strict schema enforcement
    to ensure compatibility with tool-generated outputs that may
    change across versions.

    Parameters
    ----------
    path : pathlib.Path
        Path to the input file.
    delimiter : str, optional
        Field delimiter (default: tab).
    comment : str, optional
        Comment prefix used to skip metadata lines.

    Returns
    -------
    list[dict[str, str]]
        List of records represented as dictionaries where
        column indices are used as keys.
    """

    records: List[Dict[str, str]] = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(comment) or not line.strip():
                continue

            parts = line.rstrip("\n").split(delimiter)

            records.append({str(i): value for i, value in enumerate(parts)})

    return records


def assemble_smorf_predictions(
    *,
    genome_dir: Path,
) -> List[Dict[str, str]]:
    """
    Assemble final smORF prediction records from smORFinder outputs.

    This function merges information from three intermediate files:

    - model_predictions.tsv
    - prodigal.small.gff
    - hmmsearch.tbl

    The resulting records contain coordinates, sequence context,
    probability scores, and optional smORFam annotations.

    Parameters
    ----------
    genome_dir : pathlib.Path
        Path to the genome-level smORFinder output directory.

    Returns
    -------
    list[dict[str, str]]
        Combined smORF prediction records.
    """

    tmp_dir = genome_dir / "tmp"

    # Load model predictions
    model_records = _parse_model_predictions(tmp_dir / "model_predictions.tsv")

    # Parse prodigal.small.gff
    gff_records = _parse_prodigal_small_gff(tmp_dir / "prodigal.small.gff")

    # Parse hmmsearch.tbl
    hmm_records = _parse_hmmsearch_tbl(tmp_dir / "hmmsearch.tbl")

    # Merge all sources
    gff_dict = {r["seqid"]: r for r in gff_records}
    merged: List[Dict[str, str]] = []

    for rec in model_records:
        seqid = rec.get("seqid", "")
        merged_rec = {
            "seqid": seqid,
            "contig": gff_dict.get(seqid, {}).get("contig", ""),
            "start": gff_dict.get(seqid, {}).get("start", ""),
            "end": gff_dict.get(seqid, {}).get("end", ""),
            "orient": gff_dict.get(seqid, {}).get("orient", ""),
            "smorfam": hmm_records.get(seqid, {}).get("smorfam", ""),
            "hmm_smorfam_evalue": hmm_records.get(seqid, {}).get(
                "hmm_smorfam_evalue", ""
            ),
            "dsn1_prob_smorf": rec.get("dsn1_prob_smorf", ""),
            "dsn2_prob_smorf": rec.get("dsn2_prob_smorf", ""),
            "5p_seq": rec.get("5p_seq", ""),
            "orf": rec.get("orf", ""),
            "3p_seq": rec.get("3p_seq", ""),
        }
        merged.append(merged_rec)

    return merged


# Internal Helpers
def _parse_model_predictions(
    path: Path,
) -> List[Dict[str, str]]:
    """
    Parse smORFinder model prediction table.

    Parameters
    ----------
    path : pathlib.Path
        Path to model_predictions.tsv.

    Returns
    -------
    list[dict[str, str]]
        Parsed prediction records.
    """

    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        model_records = [dict(r) for r in reader]

    for r in model_records:
        if "orf_seq" in r:
            r["orf"] = r.pop("orf_seq")  # Standardize field name

    return model_records


def _parse_prodigal_small_gff(
    path: Path,
) -> List[Dict[str, str]]:
    """
    Parse Prodigal GFF output for small ORFs.

    Parameters
    ----------
    path : pathlib.Path
        Path to prodigal.small.gff.

    Returns
    -------
    list[dict[str, str]]
        Parsed coordinate records.
    """

    gff_records: List[Dict[str, str]] = []

    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue

            seqid = ""
            attrs = cols[8]
            if "ID=" in attrs:
                seqid = attrs.split("ID=")[1].split(";")[0]

            gff_records.append(
                {
                    "seqid": seqid,
                    "contig": cols[0],
                    "start": cols[3],
                    "end": cols[4],
                    "orient": cols[6],
                }
            )

    return gff_records


def _parse_hmmsearch_tbl(
    path: Path,
) -> Dict[str, Dict[str, str]]:
    """
    Parse HMMER hmmsearch table output.

    Parameters
    ----------
    path : pathlib.Path
        Path to hmmsearch.tbl.

    Returns
    -------
    dict[str, dict[str, str]]
        Mapping from sequence ID to HMM annotation.
    """

    hmm_records: Dict[str, Dict[str, str]] = {}

    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split()
            seqid, smorfam, evalue = cols[0], cols[2], cols[4]

            # Keep the lowest e-value per sequence ID
            if seqid not in hmm_records or float(evalue) < float(
                hmm_records[seqid]["hmm_smorfam_evalue"]
            ):
                hmm_records[seqid] = {
                    "smorfam": smorfam,
                    "hmm_smorfam_evalue": evalue,
                }

    return hmm_records
