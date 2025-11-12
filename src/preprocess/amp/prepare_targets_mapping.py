# pylint: disable=import-error, wrong-import-position, too-many-locals
"""
Apply manually curated YAML mapping to AMP targets table.
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
from src.utils.io_utils import file_exists, load_dataframe_by_columns, load_yaml_mapping
from src.utils.log_utils import (
    get_pipeline_completion_message,
    get_task_completion_message,
)


# ============================== Custom Functions ==============================
def manual_targets_mapping(
    input_yaml: str,
    amp_input_csv: str,
    logger: CustomLogger,
) -> None:
    """
    Apply manually curated YAML mapping to AMP targets CSV.

    Parameters
    ----------
    input_yaml : str
        Path to YAML mapping file (manual corrections).
    amp_input_csv : str
        Path to cleaned AMP targets CSV.
    logger : CustomLogger
        Logger instance to track progress.

    Returns
    -------
    None
    """
    logger.info("/ Task: Apply manual target mappings to normalized table")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Load YAML mapping and summarize its composition
        manual_mapping = load_yaml_mapping(file_path=input_yaml)
        yaml_total = len(manual_mapping)
        yaml_unknown = sum(
            1 for v in manual_mapping.values() if str(v).lower() == "unknown"
        )
        yaml_mapped = yaml_total - yaml_unknown
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ YAML Mapping Overview ]\n"
                f"▸ Total YAML entries : {yaml_total}\n"
                f"   - 'Valid'   : {yaml_mapped}\n"
                f"   - 'Unknown' : {yaml_unknown}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Load AMP targets CSV and apply YAML-based normalization
        df = load_dataframe_by_columns(
            file_path=amp_input_csv,
            required_columns=["ID", "Sequence", "Targets"],
            has_header=True,
        )
        df["Mapped_Targets"] = df["Targets"].map(manual_mapping)

        # Summary statistics
        total_records = len(df)
        mapped_col = df["Mapped_Targets"].astype(str).str.strip().str.lower()
        unknown_mask = mapped_col.eq("unknown")
        unmapped_mask = df["Mapped_Targets"].isna()
        mapped_mask = ~unknown_mask & ~unmapped_mask
        mapped_count = mapped_mask.sum()
        unknown_count = unknown_mask.sum()
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Post-Mapping Summary ]\n"
                f"▸ Total records       : {total_records}\n"
                f"▸ Mapped successfully : {mapped_count}\n"
                f"▸ Marked as Unknown   : {unknown_count}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Save results
        df = df.drop(columns=["Targets"])
        df = df.rename(columns={"Mapped_Targets": "Targets"})
        base_dir = os.path.dirname(amp_input_csv)
        output_path = os.path.join(base_dir, "targets_normalized.csv")
        df.to_csv(output_path, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Saved:\n'{output_path}'",
            border="|",
            length=120,
        )

        # Export unmapped and unknown entries for inspection
        unmapped_df = df[unmapped_mask].copy()
        if not unmapped_df.empty:
            unmapped_output = os.path.join(base_dir, "targets_unmapped.csv")
            unmapped_df.to_csv(unmapped_output, index=False)
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    f"[ Unmapped Targets Exported ]\n" f"▸ Entries: {len(unmapped_df)}"
                ),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            logger.log_with_borders(
                level=logging.INFO,
                message=f"Saved:\n'{unmapped_output}'",
                border="|",
                length=120,
            )

        # Log completion
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error during 'apply_manual_targets_mapping()'.")
        raise


# ============================== Pipeline Entry Point ==============================
def run_normalize_targets(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Run the pipeline step for applying manually curated YAML mappings
    to the cleaned AMP targets table.

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
        input_yaml = kwargs.get(
            "mapping_input_yaml",
            os.path.join(base_path, "data/manual/targets_mapping/targets_mapping.yml"),
        )
        amp_input_csv = kwargs.get(
            "mapping_amp_input_csv",
            os.path.join(base_path, "data/processed/dbAMP/targets_clean.csv"),
        )

        # Check input files exist
        for path in [input_yaml, amp_input_csv]:
            if not file_exists(file_path=path):
                raise FileNotFoundError(f"File not found: '{path}'")

        # Execute Mapping
        manual_targets_mapping(
            input_yaml=input_yaml,
            amp_input_csv=amp_input_csv,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_normalize_targets()'")
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
