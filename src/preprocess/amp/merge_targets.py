# pylint: disable=line-too-long, import-error, wrong-import-position, too-many-locals
"""
Merge normalized AMP targets from dbAMP, DBAASP, and DRAMP by sequence.
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import sys
import time
from collections import Counter, defaultdict
from typing import DefaultDict, Set, TypedDict

# ============================== Third-Party Library Imports ==============================
import numpy as np
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
class MergeBucket(TypedDict):
    """
    Typed dictionary representing a single merged AMP entry.

    Keys
    -----
    IDs : Set[str]
        All AMP IDs associated with the same sequence across databases.
    Targets : Set[str]
        All unique target organisms or cell lines associated with the sequence.
    RepID : str
        The representative ID chosen based on priority rule (dbAMP > DBAASP > DRAMP).
    Length : int
        Peptide sequence length.
    n_Targets : int
        Number of unique targets associated with the sequence.
    """

    IDs: Set[str]
    Targets: Set[str]
    RepID: str
    Length: int
    n_Targets: int


def _bucket_factory() -> MergeBucket:
    return {
        "IDs": set(),
        "Targets": set(),
        "RepID": "",
        "Length": 0,
        "n_Targets": 0,
    }


def _id_priority(id_str: str) -> int:
    """
    Map an AMP ID to a numeric priority (lower is better).

    Parameters
    ----------
    id_str : str
        Raw AMP identifier.

    Returns
    -------
    int
        Priority bucket. Lower value indicates higher priority.
    """
    s = id_str.strip().lower()
    if s.startswith("dbamp_"):
        return 0
    if s.startswith("dbaasp_"):
        return 1
    if s.startswith("dramp_") or s.startswith("dramp"):
        return 2
    return 3


def _load_and_filter_unknown(path: str, name: str) -> pd.DataFrame:
    """
    Load a normalized AMP target CSV and remove entries with 'Unknown' or empty targets.
    Raise an error if any NaN values are detected in 'Targets'.

    Parameters
    ----------
    path : str
        Path to the input CSV file.
    name : str
        Name of the source database (for error message context).

    Returns
    -------
    pd.DataFrame
        Cleaned DataFrame with valid target entries only.
    """
    df = load_dataframe_by_columns(
        file_path=path, required_columns=["ID", "Sequence", "Targets"]
    )

    # Strict NaN check
    if df["Targets"].isna().any():
        nan_count = df["Targets"].isna().sum()
        raise ValueError(
            f"[{name}] Detected {nan_count} NaN entries in 'Targets' column. "
            "This indicates unexpected missing target data in the normalized CSV."
        )

    # Remove 'Unknown' and empty string entries
    mask_valid = (df["Targets"].astype(str).str.strip().str.lower() != "unknown") & (
        df["Targets"].astype(str).str.strip() != ""
    )
    df = df[mask_valid].copy()

    return df


def choose_rep_id(ids: Set[str]) -> str:
    """
    Choose a representative ID from a set of IDs using a stable rule.

    Parameters
    ----------
    ids : Set[str]
        Set of AMP IDs collected for a sequence.
    """
    if not ids:
        return ""
    return sorted((i.strip() for i in ids), key=lambda x: (_id_priority(x), x))[0]


def normalize_seq(seq: str) -> str:
    """
    Normalize peptide sequence into a canonical uppercase format.

    Parameters
    ----------
    seq : str
        Raw sequence string extracted from input CSV file.
        May contain lowercase letters, trailing spaces, or NaN values.

    Returns
    -------
    str
        Canonicalized peptide sequence string.
        Returns an empty string ("") if the input is NaN or invalid.
    """
    if pd.isna(seq):
        return ""
    return str(seq).strip().upper()


def log_sequence_length_stats(merged: dict[str, dict], logger: CustomLogger) -> None:
    """
    Compute and log sequence length statistics.

    Parameters
    ----------
    merged : dict
        Dictionary containing merged AMP data with sequence keys and metadata values.
    logger : CustomLogger
        Logger instance for structured logging.

    Returns
    -------
    None
    """
    all_lengths = np.array([bucket["Length"] for bucket in merged.values()])
    unique_seq_count = len(all_lengths)

    # Core descriptive stats
    min_len = int(all_lengths.min())
    max_len = int(all_lengths.max())
    mean_len = float(all_lengths.mean())
    median_len = float(np.median(all_lengths))
    std_len = float(all_lengths.std())
    q1, q3 = np.percentile(all_lengths, [25, 75])
    iqr = q3 - q1
    skewness = float(((all_lengths - mean_len) ** 3).mean() / (std_len**3))
    kurtosis = float(((all_lengths - mean_len) ** 4).mean() / (std_len**4))

    # Most frequent lengths
    length_counts = Counter(all_lengths)
    top_lengths = length_counts.most_common(5)
    top_lengths_str = ", ".join([f"{int(k)}aa ({v})" for k, v in top_lengths])

    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
    logger.log_with_borders(
        level=logging.INFO,
        message=(
            "[ Sequence Length Statistics ]\n"
            f"▸ Total unique sequences             : {unique_seq_count}\n"
            f"▸ Length (min / mean / median / max) : {min_len} / {mean_len:.2f} / {median_len:.2f} / {max_len}\n"
            f"▸ Standard deviation (SD)            : {std_len:.2f}\n"
            f"▸ 25th / 75th percentile             : {q1:.1f} / {q3:.1f}\n"
            f"▸ Interquartile range (IQR)          : {iqr:.1f}\n"
            f"▸ Skewness / Kurtosis                : {skewness:.2f} / {kurtosis:.2f}\n"
            f"▸ Most frequent lengths              : {top_lengths_str}"
        ),
        border="|",
        length=120,
    )


def log_target_stats(merged: dict[str, dict], logger: CustomLogger) -> None:
    """
    Compute and log target distribution statistics.

    Parameters
    ----------
    merged : dict
        Dictionary containing merged AMP data with sequence keys and metadata values.
    logger : CustomLogger
        Logger instance for structured logging.

    Returns
    -------
    None
    """
    all_target_lists = [bucket["Targets"] for bucket in merged.values()]
    target_counts = np.array([len(tset) for tset in all_target_lists])
    unique_targets = set().union(*all_target_lists)

    # Basic stats
    unique_target_count = len(unique_targets)
    avg_targets = float(target_counts.mean())
    median_targets = float(np.median(target_counts))
    max_targets = int(target_counts.max())
    multi_target = int((target_counts > 1).sum())
    high_multi_target = int((target_counts >= 10).sum())

    # Most frequent targets
    target_counter = Counter()
    for tset in all_target_lists:
        target_counter.update(tset)
    top_targets = target_counter.most_common(10)
    top_targets_str = "\n  - ".join([f"'{t}' ({c})" for t, c in top_targets])

    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
    logger.log_with_borders(
        level=logging.INFO,
        message=(
            "[ Target Statistics ]\n"
            f"▸ Unique targets             : {unique_target_count}\n"
            f"▸ Avg targets per sequence   : {avg_targets:.2f}\n"
            f"▸ Median targets per seq     : {median_targets:.1f}\n"
            f"▸ Max targets per sequence   : {max_targets}\n"
            f"▸ Multi-target sequences     : {multi_target}\n"
            f"▸ Sequences with ≥10 targets : {high_multi_target}\n"
            f"▸ Most frequent targets      :\n  - {top_targets_str}"
        ),
        border="|",
        length=120,
    )


def merge_targets_by_sequence(
    dbamp_csv: str,
    dbaasp_csv: str,
    dramp_csv: str,
    output_dir: str,
    logger: CustomLogger,
) -> None:
    """
    Merge AMP target data from dbAMP, DBAASP, and DRAMP by unique Sequence.

    Parameters
    ----------
    dbamp_csv : str
        Path to dbAMP targets CSV file.
    dbaasp_csv : str
        Path to DBAASP targets CSV file.
    dramp_csv : str
        Path to DRAMP targets CSV file.
    output_dir : str
        Directory to save the merged CSV and FASTA outputs.
    logger : CustomLogger
        Logger instance for structured logging.

    Returns
    -------
    None
    """
    logger.info("/ Task: Merge AMP targets from three databases by 'Sequence'")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Load input data
        df_dbamp = _load_and_filter_unknown(dbamp_csv, "dbAMP")
        df_dbaasp = _load_and_filter_unknown(dbaasp_csv, "DBAASP")
        df_dramp = _load_and_filter_unknown(dramp_csv, "Dramp")
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Input Summary ]\n"
                f"▸ dbAMP  : {len(df_dbamp)} records\n"
                f"▸ DBAASP : {len(df_dbaasp)} records\n"
                f"▸ DRAMP  : {len(df_dramp)} records"
            ),
            border="|",
            length=120,
        )

        # Initialize Merge Dictionary
        merged: DefaultDict[str, MergeBucket] = defaultdict(_bucket_factory)

        # Populate Dictionary
        for df in [df_dbamp, df_dbaasp, df_dramp]:
            for _, row in df.iterrows():
                seq = normalize_seq(row["Sequence"])
                if not seq:
                    continue
                id_str = str(row["ID"]).strip()
                tgt = str(row["Targets"]).strip()

                # collect
                bucket = merged[seq]
                bucket["IDs"].add(id_str)
                bucket["Targets"].add(tgt)

        # Post-process: compute RepID per sequence using the priority rule
        for seq, bucket in merged.items():
            bucket["RepID"] = choose_rep_id(bucket["IDs"])
            bucket["Length"] = len(seq)
            bucket["n_Targets"] = len(bucket["Targets"])

        # Log statistical summaries
        log_sequence_length_stats(merged, logger)  # type: ignore
        log_target_stats(merged, logger)  # type: ignore

        # Save Outputs (long-format CSV + FASTA)
        records = []
        for seq, bucket in merged.items():
            rep_id = bucket["RepID"]
            for target in sorted(bucket["Targets"]):
                records.append(
                    {
                        "ID": rep_id,
                        "Sequence": seq,
                        "Targets": target,
                    }
                )

        df_out = pd.DataFrame(records)

        # Ensure output directory
        os.makedirs(output_dir, exist_ok=True)

        # Define paths
        csv_long_path = os.path.join(output_dir, "targets_merged.csv")
        csv_unique_path = os.path.join(output_dir, "targets_unique.csv")
        fasta_path = os.path.join(output_dir, "targets_merged.fasta")

        # Save long-format CSV
        df_out.to_csv(csv_long_path, index=False)

        # Save FASTA (each RepID once)
        with open(fasta_path, "w", encoding="utf-8") as f:
            for rep_id, seq in df_out.drop_duplicates("ID")[
                ["ID", "Sequence"]
            ].itertuples(index=False):
                f.write(f">{rep_id}\n{seq}\n")

        # Save unique Targets list
        unique_targets = sorted(df_out["Targets"].unique())
        df_unique = pd.DataFrame(unique_targets, columns=["Targets"])
        df_unique.to_csv(csv_unique_path, index=False)

        # Log summary
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        logger.log_with_borders(
            level=logging.INFO,
            message=(f"Saved:\n'{csv_long_path}'\n'{csv_unique_path}'\n'{fasta_path}'"),
            border="|",
            length=120,
        )

        # Completion
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error during 'merge_targets_by_sequence()'.")
        raise


# ============================== Pipeline Entry ==============================
def run_merge_targets_by_sequence(
    base_path: str, logger: CustomLogger, **kwargs
) -> None:
    """
    Merge AMP targets by sequence across all databases.

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
        dbamp_csv = kwargs.get(
            "merge_dbamp_csv",
            os.path.join(base_path, "data/processed/dbAMP/targets_normalized.csv"),
        )
        dbaasp_csv = kwargs.get(
            "merge_dbaasp_csv",
            os.path.join(base_path, "data/processed/DBAASP/targets_normalized.csv"),
        )
        dramp_csv = kwargs.get(
            "merge_dramp_csv",
            os.path.join(base_path, "data/processed/DRAMP/targets_normalized.csv"),
        )
        output_dir = kwargs.get(
            "merge_output_dir",
            os.path.join(base_path, "data/processed/merged/"),
        )

        # Validate input file
        for path in [dbamp_csv, dbaasp_csv, dramp_csv]:
            if not file_exists(path):
                raise FileNotFoundError(f"File not found: {path}")

        # Execute
        merge_targets_by_sequence(
            dbamp_csv=dbamp_csv,
            dbaasp_csv=dbaasp_csv,
            dramp_csv=dramp_csv,
            output_dir=output_dir,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_merge_targets_by_sequence()'")
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
