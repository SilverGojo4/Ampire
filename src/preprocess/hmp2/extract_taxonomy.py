# pylint: disable=import-error, wrong-import-position,
"""
Extract genus-level taxonomy from QIIME2 (16S) OTU/ASV abundance tables.
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
from src.utils.io_utils import directory_exists, file_exists
from src.utils.log_utils import (
    get_pipeline_completion_message,
    get_task_completion_message,
)


# ============================== Custom Functions ==============================
def extract_genus_from_string(tax_string: str) -> str | None:
    """
    Given a taxonomy string ('A;B;C;D;E;F'), extract the genus-level (6th layer).

    Parameters
    ----------
    tax_string : str
        Full taxonomy string.

    Returns
    -------
    str | None
        Genus name if taxonomy reaches depth >= 6; otherwise None.
    """
    if not isinstance(tax_string, str):
        return None

    parts = [p for p in tax_string.split(";") if p]
    if len(parts) < 6:
        return None

    return parts[5].strip()


def sort_genus_list(genus_list: list[str]) -> list[str]:
    """
    Custom sorting:
    1. Normal genus A–Z first
    2. Bracketed genus (e.g., '[Eubacterium] eligens group') placed last

    Parameters
    ----------
    genus_list : list[str]

    Returns
    -------
    list[str]
    """

    def custom_key(name):
        return (1 if name.startswith("[") else 0, name.lower())

    return sorted(genus_list, key=custom_key)


def extract_hmp2_genus_taxonomy(
    input_tsv: str,
    output_csv: str,
    logger: CustomLogger,
) -> None:
    """
    Extract genus from QIIME2 output abundance table (TSV) and save sorted unique genus list.

    Parameters
    ----------
    input_tsv : str
        Path to QIIME2 relative abundance .tsv (should contain '#OTU ID').
    output_csv : str
        Path to save unique genus list as CSV.
    logger : CustomLogger
        Logger instance.

    Returns
    -------
    None
    """
    logger.info("/ Task: Extract genus-level taxonomy from HMP2 OTU tables")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Load relevant columns
        df = pd.read_csv(input_tsv, sep="\t", skiprows=1)
        if "#OTU ID" not in df.columns:
            raise ValueError("Column '#OTU ID' not found in input TSV.")

        # Extract genus
        taxonomy_raw = df["#OTU ID"].astype(str)
        genus_series = taxonomy_raw.apply(extract_genus_from_string)
        genus_clean = genus_series.dropna().unique().tolist()
        genus_sorted = sort_genus_list(genus_clean)
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Genus Extraction Summary ]\n"
                f"▸ Extracted entries  : {len(genus_series)}\n"
                f"▸ Unique genus count : {len(genus_sorted)}\n"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Write output
        output_dir = os.path.dirname(output_csv)
        if not directory_exists(dir_path=output_dir):
            os.makedirs(output_dir)
        pd.DataFrame({"Genus": genus_sorted}).to_csv(output_csv, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Saved:\n'{output_csv}'",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Summary
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error in 'extract_hmp2_genus_taxonomy()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_extract_hmp2_genus_taxonomy(
    base_path: str, logger: CustomLogger, **kwargs
) -> None:
    """
    Pipeline entry point for extracting genus-level taxonomy from OTU tables.

    Parameters
    ----------
    base_path : str
        Project root path.
    logger : CustomLogger
        Logger instance.
    **kwargs : dict
        Optional CLI overrides.

    Returns
    -------
    None
    """
    # Start timing
    start_time = time.time()

    try:
        # -------------------- Retrieve input parameters (CLI or default) --------------------
        input_tsv = kwargs.get(
            "hmp2_rel_abundance_tsv",
            os.path.join(
                base_path,
                "experiments/HMP2/biopsy_16S/qiime2/abundance_tables/abs-abund-table-6.tsv",
            ),
        )
        output_csv = kwargs.get(
            "hmp2_genus_output_csv",
            os.path.join(base_path, "data/processed/HMP2/genus_list.csv"),
        )

        # Check input file exists
        if not file_exists(file_path=input_tsv):
            raise FileNotFoundError(f"File not found: '{input_tsv}'")

        # Execute extraction
        extract_hmp2_genus_taxonomy(
            input_tsv=input_tsv,
            output_csv=output_csv,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_extract_hmp2_genus_taxonomy()'")
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
