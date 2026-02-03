# pylint: disable=
"""
Genus-level input loaders for Ampire pipelines.

This module defines utilities for parsing and normalizing
genus-level input specifications used by batch pipelines.
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import List

# Third-Party Library Imports
import pandas as pd


# Public API
def load_genus_list(
    genus_list_path: Path,
) -> List[str]:
    """
    Load and normalize a list of genus names for batch pipeline execution.

    This function parses human-defined input specifications and converts
    them into a deterministic, pipeline-ready list of genus names.

    Supported formats
    ------------------
    - .txt : one genus name per line
             lines starting with '#' are treated as comments
    - .csv : must contain a column named 'genus'

    Parameters
    ----------
    genus_list_path : pathlib.Path
        Path to the genus list file.

    Returns
    -------
    list[str]
        Normalized and deduplicated genus names, preserving
        original ordering.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If the file format is unsupported or required
        columns are missing.
    """
    if not genus_list_path.exists():
        raise FileNotFoundError(f"Genus list not found: {genus_list_path}")

    suffix = genus_list_path.suffix.lower()

    # TXT format: one genus per line
    if suffix == ".txt":
        genera = [
            line.strip()
            for line in genus_list_path.read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
        return _normalize_genus_list(genera)

    # CSV format: must contain 'genus' column
    if suffix == ".csv":
        df = pd.read_csv(genus_list_path)

        if "genus" not in df.columns:
            raise ValueError("CSV genus list must contain a column named 'genus'.")

        genera = df["genus"].dropna().astype(str).tolist()
        return _normalize_genus_list(genera)

    raise ValueError(
        f"Unsupported genus list format: '{suffix}'. " "Supported formats: .txt, .csv"
    )


# Internal Helpers
def _normalize_genus_list(
    genera: List[str],
) -> List[str]:
    """
    Normalize and deduplicate genus names.

    Normalization rules
    -------------------
    - strip leading/trailing whitespace
    - remove empty entries
    - preserve original capitalization
    - deduplicate while preserving input order

    Parameters
    ----------
    genera : list[str]
        Raw genus name list.

    Returns
    -------
    list[str]
        Cleaned and unique genus names.
    """
    seen = set()
    normalized: List[str] = []

    for genus in genera:
        name = genus.strip()
        if not name:
            continue
        if name in seen:
            continue

        seen.add(name)
        normalized.append(name)

    return normalized
