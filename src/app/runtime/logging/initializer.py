# pylint: disable=
"""
Runtime logging initialization utilities.

This module provides a standardized entry point for initializing
logging behavior for a single Ampire pipeline execution.
"""

from __future__ import annotations

# Standard Library Imports
from datetime import datetime
from pathlib import Path
from typing import Optional

# External Logging Toolkit
from logging_toolkit.logger import CustomLogger
from logging_toolkit.setup_logging import setup_logging

# App Imports
from app.utils.config_utils import load_json


# Public API
def initialize_pipeline_logger(
    *,
    pipeline_name: str,
    logging_config_path: Path,
    logs_root: Path = Path("logs"),
    run_id: Optional[str] = None,
) -> CustomLogger:
    """
    Initialize a pipeline-scoped logger for a single execution run.

    This function is the *only* place where logging configuration
    is applied during an Ampire execution.

    The resulting log directory structure follows:

        logs/{pipeline_name}/{run_id}.log

    Parameters
    ----------
    pipeline_name : str
        Pipeline identifier (e.g. "genomes.genus").
    logging_config_path : pathlib.Path
        Path to JSON logging configuration file.
    logs_root : pathlib.Path, optional
        Root directory for runtime logs.
    run_id : str, optional
        Explicit execution identifier. If omitted, a timestamp-based
        identifier will be generated automatically.

    Returns
    -------
    CustomLogger
        Logger instance bound to the current pipeline execution.
    """

    # Step 1. Normalize pipeline name
    pipeline_dirname = _normalize_pipeline_name(pipeline_name)

    # Step 2. Resolve run identifier
    if run_id is None:
        run_id = _generate_run_id()

    # Step 3. Construct log directory
    log_dir = logs_root / pipeline_dirname
    log_dir.mkdir(parents=True, exist_ok=True)

    # Step 4. Resolve final log file path
    log_file_path = log_dir / f"{run_id}.log"

    # Step 5. Validate logging config exists
    load_json(logging_config_path)

    # Step 6. Initialize logging backend
    logger = setup_logging(
        input_config_file=str(logging_config_path),
        logger_name=pipeline_name,
        handler_name="file",
        output_log_path=str(log_file_path),
    )

    return logger


# Internal Helpers
def _generate_run_id() -> str:
    """
    Generate a timestamp-based execution identifier.

    Format
    ------
    YYYYMMDD-HHMMSS

    Returns
    -------
    str
        Runtime-safe identifier for log file naming.
    """
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def _normalize_pipeline_name(pipeline_name: str) -> str:
    """
    Normalize pipeline identifiers for filesystem compatibility.

    Examples
    --------
    genomes.genus   -> genomes.genus
    genomes/genus   -> genomes.genus
    genomes\\genus  -> genomes.genus

    Parameters
    ----------
    pipeline_name : str
        Raw pipeline identifier.

    Returns
    -------
    str
        Normalized pipeline directory name.
    """
    return pipeline_name.replace("/", ".").replace("\\", ".").strip()
