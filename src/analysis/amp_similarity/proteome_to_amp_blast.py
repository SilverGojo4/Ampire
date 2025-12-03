# pylint: disable=import-error, wrong-import-position, too-many-arguments, too-many-positional-arguments, too-many-locals, broad-exception-caught
"""
BLAST proteome-derived peptides against AMP database (genus-wide pipeline)
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import sys
import time
from typing import List, Tuple

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
from src.analysis.amp_similarity.blastp_runner import blastp_proteome_to_amp
from src.utils.log_utils import (
    get_pipeline_completion_message,
    get_task_completion_message,
)


# ============================== Custom Functions ==============================
def verify_proteome_data_availability(
    genus_name: str,
    proteomes_root: str,
    logger: CustomLogger,
) -> Tuple[List[str], List[str], bool, bool, bool]:
    """
    Verify the availability of UniProt proteome FASTA files (reviewed + unreviewed)
    for a given genus.

    Parameters
    ----------
    genus_name : str
        Target genus name (e.g., "Abiotrophia").
    proteomes_root : str
        Root directory containing all genus-level proteome folders.
    logger : CustomLogger
        Logger instance for structured tracking.

    Returns
    -------
    (reviewed_fastas, unreviewed_fastas, has_reviewed, has_unreviewed, success)
        reviewed_fastas   : list of absolute paths to reviewed FASTA files
        unreviewed_fastas : list of absolute paths to unreviewed FASTA files
        has_reviewed      : True if reviewed FASTA files exist
        has_unreviewed    : True if unreviewed FASTA files exist
        success           : True if at least one of the two categories exists
    """
    logger.info(f"/ Task: Validate presence of proteome FASTA data for '{genus_name}'")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Directory Setup
        genus_dir = os.path.join(proteomes_root, genus_name)
        fasta_dir = os.path.join(genus_dir, "fasta")
        reviewed_dir = os.path.join(fasta_dir, "reviewed")
        unreviewed_dir = os.path.join(fasta_dir, "unreviewed")

        reviewed_fastas = []
        unreviewed_fastas = []
        has_reviewed = False
        has_unreviewed = False

        # Check genus folder existence
        if not os.path.exists(genus_dir):
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    f"[ Proteome data summary ]\n"
                    f"▸ Genus  : '{genus_name}'\n"
                    f"▸ Reason : No proteome directory found for this genus"
                ),
                border="|",
                length=120,
            )
        else:
            # Reviewed FASTA scanning
            if os.path.exists(reviewed_dir):
                reviewed_fastas = sorted(
                    [
                        os.path.join(reviewed_dir, f)
                        for f in os.listdir(reviewed_dir)
                        if f.endswith(".fasta")
                    ]
                )
                has_reviewed = len(reviewed_fastas) > 0

            # Unreviewed FASTA scanning
            if os.path.exists(unreviewed_dir):
                unreviewed_fastas = sorted(
                    [
                        os.path.join(unreviewed_dir, f)
                        for f in os.listdir(unreviewed_dir)
                        if f.endswith(".fasta")
                    ]
                )
                has_unreviewed = len(unreviewed_fastas) > 0

            # Summary Logging
            if not has_reviewed and not has_unreviewed:
                logger.log_with_borders(
                    level=logging.INFO,
                    message=(
                        f"[ Proteome data summary ]\n"
                        f"▸ Genus  : '{genus_name}'\n"
                        f"▸ Reason : No reviewed or unreviewed FASTA files found."
                    ),
                    border="|",
                    length=120,
                )
            else:
                logger.log_with_borders(
                    level=logging.INFO,
                    message=(
                        f"[ Proteome data summary ]\n"
                        f"▸ Genus             : '{genus_name}'\n"
                        f"▸ reviewed FASTA    : {len(reviewed_fastas)} file(s)\n"
                        f"▸ unreviewed FASTA  : {len(unreviewed_fastas)} file(s)\n"
                        f"▸ Status            : Proteome data available."
                    ),
                    border="|",
                    length=120,
                )

        # Determine overall success
        success = has_reviewed or has_unreviewed

        # Log completion
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        return reviewed_fastas, unreviewed_fastas, has_reviewed, has_unreviewed, success

    except Exception:
        logger.exception("Unexpected error in 'verify_proteome_data_availability()'")
        raise


def blast_proteome_single_batch(
    fasta_files: List[str],
    genus_name: str,
    category: str,  # "reviewed" or "unreviewed"
    reference_db: str,
    logger: CustomLogger,
    threads: int = 8,
    evalue: float = 1e-2,
    word_size: int = 3,
) -> None:
    """
    Run BLASTP for UniProt proteome FASTA files (reviewed or unreviewed)
    under a given genus.

    Parameters
    ----------
    fasta_files : List[str]
        List of FASTA file paths (reviewed or unreviewed).
    genus_name : str
        Target genus name (e.g., "Abiotrophia").
    category : str
        Either "reviewed" or "unreviewed".
    reference_db : str
        Path prefix of AMP BLAST database.
    logger : CustomLogger
        Logger instance for structured logging.
    threads : int
        Number of CPU threads.
    evalue : float
        BLASTP e-value cutoff (default=1e-2).
    word_size : int
        BLASTP word size (default=3 for regular proteins).

    Returns
    -------
    None
    """
    logger.info(f"/ Task: Run BLASTP proteome ('{category}') to AMP for '{genus_name}'")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Summary of settings
        total = len(fasta_files)
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"[ Proteome BLASTP Screening Summary ]\n"
                f"▸ Genus            : '{genus_name}'\n"
                f"▸ Category         : '{category}'\n"
                f"▸ FASTA count      : {total}\n"
                f"▸ Reference DB     : '{reference_db}'\n"
                f"▸ Threads          : {threads}\n"
                f"▸ E-value cutoff   : {evalue}\n"
                f"▸ Word size        : {word_size}\n"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Iterate over each genome FASTA
        success_count = 0
        failed = []

        # Prepare genus-level blast root
        genus_dir = os.path.abspath(os.path.join(fasta_files[0], "../../.."))
        blast_root = os.path.join(genus_dir, "blast/proteome", category)
        os.makedirs(blast_root, exist_ok=True)

        # Iterate over FASTA files
        for fasta_path in fasta_files:
            fasta_name = os.path.splitext(os.path.basename(fasta_path))[0]

            # Folder: blast/proteome/<category>/<basename>/
            output_dir = os.path.join(blast_root, fasta_name)
            os.makedirs(output_dir, exist_ok=True)
            output_csv = os.path.join(output_dir, "amp_hits.csv")

            # Verify FASTA existence
            if not os.path.exists(fasta_path):
                failed.append(fasta_name)
                continue

            try:
                blastp_proteome_to_amp(
                    query_fasta=fasta_path,
                    reference_db=reference_db,
                    output_csv=output_csv,
                    logger=logger,
                    threads=threads,
                    evalue=evalue,
                    word_size=word_size,
                )

                # Check if output created
                if not os.path.exists(output_csv):
                    failed.append(fasta_name)
                    continue

                success_count += 1

            except Exception:
                failed.append(fasta_name)
                continue

        # Log completion summary
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"[ Proteome BLASTP Screening Completed ]\n"
                f"▸ Total FASTA files : {total}\n"
                f"▸ Successfully done : {success_count}\n"
                f"▸ Failed            : {len(failed)}\n"
                + (f"▸ Failed list       : {', '.join(failed)}" if failed else "")
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Log completion
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error in 'blast_proteome_single_batch()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_proteome_to_amp_by_genus(
    base_path: str, logger: CustomLogger, **kwargs
) -> None:
    """
    Run proteome-to-AMP BLAST screening pipeline for all proteomes under a given genus.

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
        genus_name = kwargs.get("proteomes_genus_name", "Escherichia")
        output_dir = kwargs.get(
            "proteomes_output_dir",
            os.path.join(base_path, "Escherichia"),
        )
        reference_db = kwargs.get(
            "blast_reference_db",
            os.path.join(base_path, "data/processed/AMP/blast_db/amp_db"),
        )
        threads = int(kwargs.get("blast_threads", 8))
        evalue = float(kwargs.get("blast_evalue", 1e-2))
        word_size = int(kwargs.get("blast_word_size", 3))
        proteomes_root = os.path.dirname(output_dir)

        # Ensure output directory exists
        if os.path.isdir(output_dir):
            os.makedirs(name=output_dir, exist_ok=True)

        # Verify proteome data availability
        reviewed_fastas, unreviewed_fastas, has_reviewed, has_unreviewed, success = (
            verify_proteome_data_availability(
                genus_name=genus_name,
                proteomes_root=proteomes_root,
                logger=logger,
            )
        )
        logger.add_spacer(level=logging.INFO, lines=1)

        # Run BLAST screening
        if success:
            if has_reviewed:
                blast_proteome_single_batch(
                    fasta_files=reviewed_fastas,
                    genus_name=genus_name,
                    category="reviewed",
                    reference_db=reference_db,
                    logger=logger,
                    threads=threads,
                    evalue=evalue,
                    word_size=word_size,
                )
                logger.add_spacer(level=logging.INFO, lines=1)
            if has_unreviewed:
                blast_proteome_single_batch(
                    fasta_files=unreviewed_fastas,
                    genus_name=genus_name,
                    category="unreviewed",
                    reference_db=reference_db,
                    logger=logger,
                    threads=threads,
                    evalue=evalue,
                    word_size=word_size,
                )
                logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_proteome_to_amp_by_genus()'")
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
