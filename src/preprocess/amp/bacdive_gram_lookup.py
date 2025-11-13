# pylint: disable=import-error, wrong-import-position, broad-exception-caught, too-many-locals
"""
Query BacDive API to determine Gram stain classification for AMP target organisms.
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import sys
import time

# ============================== Third-Party Library Imports ==============================
import pandas as pd
from bacdive import BacdiveClient

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
from src.utils.io_utils import file_exists, load_dataframe_by_columns, load_json_config
from src.utils.log_utils import (
    get_pipeline_completion_message,
    get_task_completion_message,
)


# ============================== Custom Functions ==============================
def extract_gram_stains(morph: dict) -> dict[str, int]:
    """
    Extract and count Gram stain information from a BacDive strain's morphology record.

    Parameters
    ----------
    morph : dict
        The 'Morphology' section of a BacDive strain record.

    Returns
    -------
    dict of str -> int
        Dictionary with counts of Gram stain results:
        {
            "positive": <count>,
            "negative": <count>
        }
        Returns {"positive": 0, "negative": 0} if no Gram stain information is present.
    """
    counts = {"positive": 0, "negative": 0}
    if not isinstance(morph, dict):
        return counts

    cell_morph = morph.get("cell morphology")

    # Case 1: single dict
    if isinstance(cell_morph, dict):
        gram = cell_morph.get("gram stain")
        if gram:
            gram_norm = str(gram).strip().lower()
            if gram_norm in counts:
                counts[gram_norm] += 1

    # Case 2: list of dicts
    elif isinstance(cell_morph, list):
        for item in cell_morph:
            if isinstance(item, dict):
                gram = item.get("gram stain")
                if gram:
                    gram_norm = str(gram).strip().lower()
                    if gram_norm in counts:
                        counts[gram_norm] += 1

    return counts


def cached_gram_lookup(
    name: str,
    client: BacdiveClient,
    gram_cache: dict[str, dict],
    logger: CustomLogger,
) -> dict:
    """
    Query BacDive for Gram stain information of a given target name,
    with caching to avoid redundant queries.

    Parameters
    ----------
    name : str
        Target name (normalized, after manual mapping).
    client : BacdiveClient
        An authenticated BacDive API client.
    gram_cache : dict
        Cache dictionary {name -> result dict}.
    logger : CustomLogger
        Logger instance for tracking progress.

    Returns
    -------
    dict
        Result dictionary with Gram stain information:
        {
            "positive": int | None,
            "negative": int | None,
            "final": str
        }

        Possible values of "final":
        - "positive"     → majority of Gram stain results are positive
        - "negative"     → majority of Gram stain results are negative
        - "unclassified" → strains found but no Gram info, or tie (0/0, equal counts)
        - "not_found"    → no strains found in BacDive
        - "error"        → API or query error occurred
    """
    # Cache hit
    if name in gram_cache:
        return gram_cache[name]

    try:
        # Query BacDive
        count = client.search(taxonomy=name)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"[Querying] '{name}' → {count} strains found",
            border="|",
            length=120,
        )

        # Case 1: No strains found
        if count == 0:
            result = {"positive": None, "negative": None, "final": "not_found"}
            gram_cache[name] = result
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            return result

        # Case 2: Iterate over retrieved strains
        positive_count, negative_count = 0, 0
        for strain in client.retrieve():
            morph = strain.get("Morphology", {})
            if isinstance(morph, dict):
                gram_counts = extract_gram_stains(morph)
                positive_count += gram_counts["positive"]
                negative_count += gram_counts["negative"]

        # Case 3: Decide final classification
        if positive_count > negative_count:
            final = "positive"
        elif negative_count > positive_count:
            final = "negative"
        else:
            final = "unclassified"

        result = {
            "positive": positive_count,
            "negative": negative_count,
            "final": final,
        }
        gram_cache[name] = result

        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"▸ Final result: '{final}'\n"
                f"  [+] Positive: {positive_count}\n"
                f"  [-] Negative: {negative_count}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        return result

    except Exception:
        # Case 4: Error during query
        logger.error(f"Error occurred while querying '{name}'")
        result = {"positive": None, "negative": None, "final": "error"}
        gram_cache[name] = result
        return result


def query_bacdive_for_targets(
    input_csv: str,
    bacdive_config: str,
    output_csv: str,
    logger: CustomLogger,
) -> None:
    """
    Query BacDive API for Gram stain information of each target organism.

    Parameters
    ----------
    input_csv : str
        Path to the CSV file containing a single column "Targets".
    bacdive_config : str
        Path to the JSON config file containing BacDive credentials (email, password).
    output_csv : str
        Path to save the query results CSV.
    logger : CustomLogger
        Logger instance for structured progress tracking.

    Returns
    -------
    None
    """
    logger.info("/ Task: Query BacDive for Gram classification")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:

        # Load Input
        df_targets = load_dataframe_by_columns(input_csv, required_columns=["Targets"])
        total_targets = len(df_targets)

        def classify_target(name: str) -> str:
            """
            Classify each target name to decide whether to query BacDive.

            Returns one of:
                - "unknown"                → explicitly Unknown
                - "homo_sapiens"           → human (not a cell line)
                - "homo_sapiens_cell_line" → human-derived cell line
                - "cell_line"              → non-human cell line
                - "queryable"              → others (likely bacterial)
            """
            name_clean = str(name).strip()
            name_lower = name_clean.lower()

            # Unknown
            if name_lower == "unknown":
                return "unknown"

            # Homo sapiens cases
            if "homo sapiens" in name_lower:
                if "cell line" in name_lower:
                    return "homo_sapiens_cell_line"
                return "homo_sapiens"

            # Other species cell lines
            if "cell line" in name_lower:
                return "cell_line"

            # Default — likely bacterial
            return "queryable"

        df_targets["skip_reason"] = df_targets["Targets"].apply(classify_target)
        df_targets["is_queryable"] = df_targets["skip_reason"] == "queryable"

        # Per-class counts (for logging)
        num_queryable = int(df_targets["is_queryable"].sum())
        num_skipped = total_targets - num_queryable
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Input Summary ]\n"
                f"▸ Total targets                 : {total_targets}\n"
                f"▸ Queryable (BacDive)           : {num_queryable}\n"
                f"▸ Skipped (unknown/human/cell)  : {num_skipped}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # BacDive Login
        config = load_json_config(bacdive_config)
        client = BacdiveClient(config["email"], config["password"])
        client.setSearchType("exact")
        logger.log_with_borders(
            level=logging.INFO,
            message="Successfully logged into BacDive API",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Query BacDive
        gram_cache: dict[str, dict] = {}
        results = []

        for _, row in df_targets.iterrows():
            name = row["Targets"]
            if not row["is_queryable"]:
                results.append(
                    {"positive": None, "negative": None, "final": row["skip_reason"]}
                )
                continue

            result = cached_gram_lookup(name, client, gram_cache, logger)
            results.append(result)

        # Build Result
        df_results = pd.DataFrame(results)
        df_final = pd.concat([df_targets, df_results], axis=1)

        # Summary Statistics
        summary_counts = df_final["final"].value_counts().to_dict()
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ BacDive Query Summary ]\n"
                + "\n".join(f"▸ {k:<12}: {v}" for k, v in summary_counts.items())
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Save Output
        df_final.to_csv(output_csv, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Saved:\n'{output_csv}'",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Completion
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error in 'query_bacdive_for_targets()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_query_bacdive_targets(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Run the pipeline step for querying BacDive to retrieve Gram stain classification
    for each target organism in the AMP dataset.

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
        input_csv = kwargs.get(
            "bacdive_input_csv",
            os.path.join(base_path, "data/processed/merged/targets_unique.csv"),
        )
        bacdive_config = kwargs.get(
            "bacdive_config",
            os.path.join(base_path, "configs/bacdive.json"),
        )
        output_csv = kwargs.get(
            "bacdive_output_csv",
            os.path.join(base_path, "data/processed/merged/targets_bacdive_gram.csv"),
        )

        # Check input files exist
        for path in [input_csv, bacdive_config]:
            if not file_exists(file_path=path):
                raise FileNotFoundError(f"File not found: '{path}'")

        # Execute BacDive query
        query_bacdive_for_targets(
            input_csv=input_csv,
            bacdive_config=bacdive_config,
            output_csv=output_csv,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_query_bacdive_targets()'")
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
