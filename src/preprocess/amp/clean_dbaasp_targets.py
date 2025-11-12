# pylint: disable=import-error, wrong-import-position
"""
Clean DBAASP 'Targets' text into a normalized long table (text cleaning only).
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import sys
import time

# ============================== Third-Party Library Imports ==============================
import pandas as pd

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
from src.utils.io_utils import directory_exists, file_exists, load_dataframe_by_columns
from src.utils.log_utils import (
    get_pipeline_completion_message,
    get_task_completion_message,
)


# ============================== Custom Functions ==============================
def filter_natural_aa(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter out sequences that do not contain only natural amino acids.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing a 'Sequence' column.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing only sequences composed of natural amino acids.
    """
    # Define the set of 20 natural amino acids
    natural_aa_set = set("ACDEFGHIKLMNPQRSTVWY")

    # Validate column existence
    if "Sequence" not in df.columns:
        raise KeyError("Input DataFrame must contain a 'Sequence' column.")

    # Add a boolean column to mark sequences with only natural amino acids
    df["is_natural"] = df["Sequence"].apply(
        lambda seq: all(char in natural_aa_set for char in str(seq))
    )

    # Keep only natural sequences
    filtered_df = df[df["is_natural"]].copy()

    # Remove helper column
    filtered_df.drop(columns=["is_natural"], inplace=True)

    return filtered_df


def clean_targets_table(
    input_path: str, output_path: str, logger: CustomLogger
) -> None:
    """
    Load DBAASP peptide data, keep only monomeric natural-AA sequences
    with valid biological targets, and export both the cleaned table
    and a unique list of targets.

    Parameters
    ----------
    input_path : str
        Path to the DBAASP CSV file.
    output_path : str
        Path to save the cleaned long-format CSV.
    logger : CustomLogger
        Logger instance to track progress.

    Returns
    -------
    None
    """
    logger.info(
        "/ Task: Normalize raw DBAASP targets into structured long-format table"
    )
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Load input data
        required_columns = [
            "ID",
            "COMPLEXITY",
            "SEQUENCE",
            "TARGET ACTIVITY - TARGET SPECIES",
        ]
        df = load_dataframe_by_columns(
            file_path=input_path, required_columns=required_columns
        )
        total = df.shape[0]

        # Retain monomeric peptides
        df = df[df["COMPLEXITY"] == "Monomer"]

        # Rename columns for consistency
        df = df.rename(
            columns={
                "SEQUENCE": "Sequence",
                "TARGET ACTIVITY - TARGET SPECIES": "Targets",
            }
        )[["ID", "Sequence", "Targets"]]

        # Standardize DBAASP_ID format with prefix
        df["ID"] = df["ID"].apply(lambda x: f"DBAASP_{int(x)}")

        # Filter natural amino acid sequences
        df = filter_natural_aa(df)
        total_filtered = df.shape[0]

        # Keep only non-empty Targets
        df = df[df["Targets"].notna()]
        df["Targets"] = df["Targets"].astype(str).str.strip()
        df = df[df["Targets"] != ""]
        kept = df.shape[0]
        dropped = total_filtered - kept

        # Ensure output directory
        outdir = os.path.dirname(output_path)
        if outdir and not directory_exists(outdir):
            os.makedirs(outdir)

        # Save outputs
        df.to_csv(output_path, index=False)
        base, ext = os.path.splitext(output_path)
        unique_path = f"{base}_unique{ext}"
        unique_df = pd.DataFrame(
            {"Targets": sorted(df["Targets"].dropna().unique())}
        ).reset_index(drop=True)
        unique_df.to_csv(unique_path, index=False)

        # Log summary
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ DBAASP Targets Cleaning Summary ]\n"
                f"▸ Input records              : {total}\n"
                f"▸ After monomer + natural AA : {total_filtered}\n"
                f"▸ Dropped (Targets NA/empty) : {dropped}\n"
                f"▸ Retained for cleaning      : {kept}\n"
                f"▸ Unique targets             : {unique_df.shape[0]}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Saved:\n'{output_path}'\n'{unique_path}'",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Completion log
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error in 'clean_targets_table()'")
        raise


# ============================== Pipeline Entry ==============================
def run_clean_targets(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Clean AMP 'Targets' from raw DBAASP table and export a normalized long CSV.

    Parameters
    ----------
    base_path : str
        Project root path.
    logger : CustomLogger
        Logger instance for structured logging.
    **kwargs : dict
        Optional overrides for input/output file paths.

    Returns
    -------
    None
    """
    # Start timing
    start_time = time.time()

    try:
        input_path = kwargs.get(
            "targets_input_path",
            os.path.join(base_path, "data/raw/DBAASP/peptides-complete1220.csv"),
        )
        output_path = kwargs.get(
            "targets_output_csv",
            os.path.join(base_path, "data/processed/DBAASP/targets_clean.csv"),
        )

        # Validate input file
        if not file_exists(input_path):
            raise FileNotFoundError(f"File not found: '{input_path}'")

        # Execute
        clean_targets_table(
            input_path=input_path,
            output_path=output_path,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_clean_targets()'")
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
