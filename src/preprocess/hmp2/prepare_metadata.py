# pylint: disable=import-error, wrong-import-position
"""
Prepare HMP2 metadata for nf-core/ampliseq
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import sys
import time

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
def extract_hmp2_16s_metadata(
    input_csv: str,
    output_tsv: str,
    logger: CustomLogger,
) -> None:
    """
    Filter HMP2 metadata for biopsy_16S samples and export standardized sample info.

    Parameters
    ----------
    input_csv : str
        Path to the input raw HMP2 metadata CSV.
    output_tsv : str
        Path to save the filtered metadata TSV.
    logger : CustomLogger
        Logger instance for tracking progress.

    Returns
    -------
    None
    """
    logger.info("/ Task: Prepare HMP2 metadata for 'nf-core/ampliseq'")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Load relevant columns
        required_columns = [
            "Project",
            "External ID",
            "Participant ID",
            "data_type",
            "biopsy_location",
        ]
        df = load_dataframe_by_columns(
            file_path=input_csv, required_columns=required_columns
        )

        # Filter for biopsy_16S
        df_filtered = df[df["data_type"] == "biopsy_16S"]

        # Keep only required columns
        df_final = df_filtered[["External ID", "biopsy_location"]].copy()

        # Add prefix 's' to External ID values
        df_final["External ID"] = df_final["External ID"].apply(lambda x: f"s{x}")

        # Rename columns
        df_final.rename(
            columns={
                "External ID": "sample_name",
                "biopsy_location": "sample_index",
            },
            inplace=True,
        )

        # Statistics logging
        total_samples = df_final.shape[0]
        unique_sites = df_final["sample_index"].nunique()
        counts_by_site = df_final["sample_index"].value_counts().to_dict()
        max_site_len = max(
            (len(str(site)) for site in counts_by_site.keys()), default=10
        )
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Metadata summary ]\n"
                f"▸ Total biopsy_16S samples : {total_samples}\n"
                f"▸ Unique biopsy sites : {unique_sites}\n"
                + "\n".join(
                    f"  - {site:<{max_site_len}s} : {count}"
                    for site, count in counts_by_site.items()
                )
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Ensure output directory exists
        output_dir = os.path.dirname(output_tsv)
        if not directory_exists(dir_path=output_dir):
            os.makedirs(name=output_dir)

        # Save to TSV
        df_final.to_csv(output_tsv, sep="\t", index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Saved:\n'{output_tsv}'",
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
        logger.exception("Unexpected error in 'extract_hmp2_16s_metadata()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_extract_hmp2_16s_metadata(
    base_path: str, logger: CustomLogger, **kwargs
) -> None:
    """
    Prepare HMP2 metadata for nf-core/ampliseq.
    Filters biopsy_16S samples and outputs standardized metadata TSV.

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
        metadata_input_csv = kwargs.get(
            "hmp2_metadata_input_csv",
            os.path.join(base_path, "data/raw/HMP2/hmp2_metadata_2018-08-20.csv"),
        )
        ampliseq_metadata_tsv = kwargs.get(
            "hmp2_metadata_output_tsv",
            os.path.join(base_path, "data/processed/HMP2/nfcore_metadata.tsv"),
        )

        # Check input file exists
        if not file_exists(file_path=metadata_input_csv):
            raise FileNotFoundError(f"File not found: '{metadata_input_csv}'")

        # Execute preparing hmp2 metadata
        extract_hmp2_16s_metadata(
            input_csv=metadata_input_csv,
            output_tsv=ampliseq_metadata_tsv,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_extract_hmp2_16s_metadata()'")
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
