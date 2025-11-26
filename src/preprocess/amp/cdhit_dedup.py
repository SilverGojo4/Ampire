# pylint: disable=import-error, wrong-import-position, too-many-arguments, too-many-positional-arguments
"""
Run CD-HIT to deduplicate AMP sequences and parse cluster membership + summary tables.
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import sys
import time
import subprocess

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
from src.utils.io_utils import file_exists, load_dataframe_by_columns
from src.utils.log_utils import (
    get_pipeline_completion_message,
    get_task_completion_message,
)


# ============================== Custom Functions ==============================
def cdhit(
    input_fasta: str,
    output_dir: str,
    identity: float,
    word_size: int,
    min_length: int,
    threads: int,
    memory_mb: int,
    logger: CustomLogger,
) -> None:
    """
    Execute CD-HIT with specified parameters.

    Parameters
    ----------
    input_fasta : str
        Input FASTA path
    output_dir : str
        Directory to store CD-HIT outputs
    identity : float
        Identity threshold (e.g., 1.0 for NR100)
    word_size : int
        CD-HIT word size (-n), depends on identity
    min_length : int
        Minimum sequence length to keep (-l). Sequences shorter than this are discarded.
    threads : int
        CD-HIT threads (-T)
    memory_mb : int
        Max memory (-M)
    logger : CustomLogger
        Logger instance for tracking.

    Returns
    -------
    None
    """
    logger.info("/ Task: Run CD-HIT for deduplication")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Ensure output directory
        os.makedirs(output_dir, exist_ok=True)

        # Define paths
        identity_tag = int(identity * 100)
        output_fasta = os.path.join(output_dir, f"amp_nr{identity_tag}.fasta")

        # Log parameters
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "Running CD-HIT with the following parameters:\n"
                f"- Input FASTA  : '{input_fasta}'\n"
                f"- Output FASTA : '{output_fasta}'\n"
                f"- Identity     : {identity}\n"
                f"- Word size    : {word_size}\n"
                f"- Min length   : {min_length}\n"
                f"- Threads      : {threads}\n"
                f"- Memory (MB)  : {memory_mb}\n"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Run CD-HIT
        cmd = [
            "cd-hit",
            "-i",
            input_fasta,
            "-o",
            output_fasta,
            "-c",
            str(identity),
            "-n",
            "5",
            "-d",
            "200",
            "-T",
            str(threads),
            "-M",
            str(memory_mb),
            "-l",
            str(min_length),
        ]
        subprocess.run(cmd, check=True)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Saved:\n'{output_dir}'",
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
        logger.exception("Unexpected error during 'cdhit()'")
        raise


def parse_cdhit_clstr(clstr_path: str) -> pd.DataFrame:
    """
    Parse a CD-HIT .clstr file and return a DataFrame with:
        Cluster_ID, ID, Is_Representative

    Parameters
    ----------
    clstr_path : str
        Full path to the .clstr file

    Returns
    -------
    pd.DataFrame
    """
    if not file_exists(clstr_path):
        raise FileNotFoundError(f".clstr file not found: {clstr_path}")

    rows = []
    current_cluster = None

    with open(clstr_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">Cluster"):
                cluster_number = line.split()[1]
                current_cluster = f"{cluster_number}"
                continue

            if not line:
                continue

            # parse: "0   60aa, >ID... *"
            try:
                _, rest = line.split(",", 1)
            except ValueError:
                continue

            rest = rest.strip()

            # Extract sequence ID
            seq_id = rest.split(">")[1].split("...")[0].strip()

            # Representative flag
            is_rep = rest.endswith("*")

            rows.append(
                {
                    "Cluster_ID": current_cluster,
                    "ID": seq_id,
                    "Is_Representative": is_rep,
                }
            )

    return pd.DataFrame(rows)


def log_cluster_stats(df_cluster: pd.DataFrame, logger: CustomLogger) -> None:
    """
    Compute and log CD-HIT cluster statistics.

    Parameters
    ----------
    df_cluster : pandas.DataFrame
        DataFrame parsed from CD-HIT .clstr file.
        Must contain columns: ["Cluster_ID", "ID"].

    logger : CustomLogger
        Logger instance for structured styled output.

    Returns
    -------
    None
        Logs summary statistics but does not return data.
    """
    # Count cluster sizes
    cluster_sizes = (
        df_cluster.groupby("Cluster_ID")["ID"].count().sort_values(ascending=False)
    )

    total_clusters = len(cluster_sizes)
    singleton_clusters = int((cluster_sizes == 1).sum())
    largest_cluster_size = int(cluster_sizes.max())

    # Top 10 largest clusters
    top10 = cluster_sizes.head(10).items()
    top10_str = "\n  - ".join([f"'Cluster {cid}': {size}" for cid, size in top10])
    logger.log_with_borders(
        level=logging.INFO,
        message=(
            "[ CD-HIT Cluster Statistics ]\n"
            f"▸ Total clusters          : {total_clusters}\n"
            f"▸ Singleton clusters      : {singleton_clusters}\n"
            f"▸ Largest cluster size    : {largest_cluster_size}\n"
            f"▸ Top 10 largest families :\n  - {top10_str}"
        ),
        border="|",
        length=120,
    )


def merge_metadata_with_clusters(
    metadata_csv: str,
    output_dir: str,
    identity: float,
    logger: CustomLogger,
) -> None:
    """
    Merge AMP metadata with CD-HIT cluster assignments and save output.

    Parameters
    ----------
    metadata_csv : str
        Path to the original AMP metadata CSV (must contain column 'ID')
    output_dir : str
        Directory where merged result will be saved
    identity : float
        Identity threshold (affects output filename)
    logger : CustomLogger
        Logger instance

    Returns
    -------
    None
    """
    logger.info("/ Task: Merge metadata with CD-HIT clustering")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:

        # Load metadata
        df_meta = load_dataframe_by_columns(file_path=metadata_csv)
        if "ID" not in df_meta.columns:
            raise ValueError("Metadata file must contain column 'ID'.")

        # Output path
        identity_tag = int(identity * 100)
        clstr_path = os.path.join(output_dir, f"amp_nr{identity_tag}.fasta.clstr")
        df_cluster = parse_cdhit_clstr(clstr_path=clstr_path)

        # Log cluster size distribution
        log_cluster_stats(df_cluster=df_cluster, logger=logger)

        # Merge into a single metadata table
        df_final = df_meta.merge(df_cluster, how="left", on="ID")

        # Define paths
        output_csv = os.path.join(output_dir, f"amp_nr{identity_tag}_metadata.csv")

        # Save Output
        df_final.to_csv(output_csv, index=False)
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
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
        logger.exception("Unexpected error during 'merge_metadata_with_clusters()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_cdhit_dedup(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Complete pipeline:

    1. Run CD-HIT
    2. Parse .clstr
    3. Merge metadata
    4. Produce membership + summary tables

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
        input_fasta = kwargs.get(
            "cdhit_input_fasta",
            os.path.join(base_path, "data/processed/AMP/raw_merged/amp_raw.fasta"),
        )
        metadata_csv = kwargs.get("cdhit_metadata_csv", None)
        output_dir = kwargs.get(
            "cdhit_output_dir",
            os.path.join(base_path, "data/processed/AMP/cdhit_nr100"),
        )
        identity = float(kwargs.get("identity", 1.0))
        word_size = int(kwargs.get("word_size", 5))
        min_length = int(kwargs.get("min_length", 5))
        threads = int(kwargs.get("threads", 8))
        memory_mb = int(kwargs.get("memory_mb", 0))

        # Check input files exist
        if not file_exists(file_path=input_fasta):
            raise FileNotFoundError(f"File not found: '{input_fasta}'")

        # Execute cdhit
        cdhit(
            input_fasta=input_fasta,
            output_dir=output_dir,
            identity=identity,
            word_size=word_size,
            min_length=min_length,
            threads=threads,
            memory_mb=memory_mb,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

        # Optional metadata merge step
        if metadata_csv and file_exists(file_path=metadata_csv):
            merge_metadata_with_clusters(
                metadata_csv=metadata_csv,
                output_dir=output_dir,
                identity=identity,
                logger=logger,
            )
            logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_cdhit_dedup()'")
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
