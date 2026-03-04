# pylint: disable=
"""
Normalize tblastn raw output into the Ampire canonical alignment table.

This module converts BLAST tabular output (outfmt 6) into a
deterministic, schema-stable alignment dataset used by downstream
Ampire genome-scale analysis stages.
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import List

# Third-Party Library Imports
import pandas as pd

# Column order emitted by tblastn (run_search outfmt 6)
RAW_COLUMNS: List[str] = [
    "qseqid",
    "qlen",
    "sseqid",
    "slen",
    "sframe",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovhsp",
    "qseq",
    "sseq",
]

# Columns requiring numeric coercion
NUMERIC_COLUMNS: List[str] = [
    "qlen",
    "slen",
    "sframe",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovhsp",
]

# Canonical Ampire alignment schema
CANONICAL_COLUMNS: List[str] = [
    "qseqid",
    "sseqid",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sframe",
    "length",
    "pident",
    "mismatch",
    "gapopen",
    "bitscore",
    "evalue",
    "qcovhsp",
    "qseq",
    "sseq",
]


# Public API
def normalize_tblastn_output(
    raw_tsv: Path,
) -> pd.DataFrame:
    """
    Normalize tblastn output into the canonical Ampire alignment table.

    Parameters
    ----------
    raw_tsv : pathlib.Path
        Path to the BLAST tabular output generated using
        outfmt 6 from tblastn searches.

    Returns
    -------
    pd.DataFrame
        Canonical alignment table with deterministic schema and
        derived alignment metrics. If the input file is empty,
        an empty table with a stable schema is returned.

    Raises
    ------
    FileNotFoundError
        If the provided tblastn output file does not exist.
    """

    raw_tsv = Path(raw_tsv)

    if not raw_tsv.exists():
        raise FileNotFoundError(f"tblastn output not found: '{raw_tsv}'")

    # Empty file is a valid discovery outcome
    if raw_tsv.stat().st_size == 0:
        return _empty_alignment_table()

    df = pd.read_csv(
        raw_tsv,
        sep="\t",
        header=None,
        names=RAW_COLUMNS,
    )

    if df.empty:
        return _empty_alignment_table()

    # Numeric coercion
    df[NUMERIC_COLUMNS] = df[NUMERIC_COLUMNS].apply(
        pd.to_numeric,
        errors="coerce",
    )

    # Deterministic stable sorting
    df = df.sort_values(
        by=[
            "qseqid",
            "sseqid",
            "sframe",
            "bitscore",
            "evalue",
            "qstart",
            "qend",
        ],
        ascending=[True, True, True, False, True, True, True],
        kind="mergesort",  # stable deterministic ordering
    ).reset_index(drop=True)

    # Canonical column ordering
    return df[CANONICAL_COLUMNS].copy()


# Internal Helpers
def _empty_alignment_table() -> pd.DataFrame:
    """
    Create an empty canonical alignment table with stable schema.

    Returns
    -------
    pd.DataFrame
        Empty alignment table preserving the canonical column
        ordering required by downstream Ampire stages.
    """
    return pd.DataFrame(columns=CANONICAL_COLUMNS)
