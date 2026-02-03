# pylint: disable=import-error
"""
Tabular Data Utility Module.

This module provides helper functions for loading and validating
tabular datasets (CSV, TSV, Excel) with explicit schema expectations.
"""

from __future__ import annotations

# Standard Library Imports
import os
from typing import List, Optional

# Third-Party Library Imports
import pandas as pd

# App Imports
from app.utils.fs_utils import is_file


# DataFrame Validation Utilities
def find_missing_columns(
    df: pd.DataFrame,
    required_columns: List[str],
) -> List[str]:
    """
    Identify required columns that are missing from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to validate.
    required_columns : list of str
        Expected column names.

    Returns
    -------
    List[str]
        Missing column names. Empty if all required columns are present.
    """
    return [col for col in required_columns if col not in df.columns]


def load_table_with_columns(
    file_path: str,
    required_columns: Optional[List[str]] = None,
    has_header: bool = True,
) -> pd.DataFrame:
    """
    Load a tabular file into a DataFrame with optional column validation.

    Supported formats
    -----------------
    - CSV (.csv)
    - TSV (.tsv)
    - Excel (.xlsx, .xls)

    Parameters
    ----------
    file_path : str
        Path to the input file.
    required_columns : list of str, optional
        Columns that must be present. If None, all columns are returned.
    has_header : bool, optional
        Whether the file includes a header row.

    Returns
    -------
    pd.DataFrame
        Loaded DataFrame.
    """
    # Check file exists
    if not is_file(file_path):
        raise FileNotFoundError(f"File not found: '{file_path}'")

    # Detect file type
    ext = os.path.splitext(file_path)[-1].lower()

    # Header mode
    header_opt = 0 if has_header else None

    # Load DataFrame
    if ext == ".csv":
        df = pd.read_csv(file_path, header=header_opt)
    elif ext == ".tsv":
        df = pd.read_csv(file_path, sep="\t", header=header_opt)
    elif ext in (".xlsx", ".xls"):
        df = pd.read_excel(file_path, header=header_opt)
    else:
        raise ValueError(f"Unsupported file extension: '{ext}'")

    # If no header, assign temporary column names
    if not has_header:
        df.columns = [f"column_{i}" for i in range(df.shape[1])]

    # Use all columns if none specified
    if required_columns is None:
        return df.copy()

    # Check required columns
    missing = find_missing_columns(df, required_columns)
    if missing:
        raise ValueError(f"Missing required columns: '{str(missing)}'")

    return df[required_columns].copy()
