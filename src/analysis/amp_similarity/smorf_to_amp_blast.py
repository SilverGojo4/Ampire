# pylint: disable=import-error, wrong-import-position, too-many-arguments, too-many-positional-arguments, too-many-locals, broad-exception-caught
"""
BLAST smORF-derived peptides against AMP database (genus-wide pipeline)
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
from src.analysis.amp_similarity.blastp_runner import blastp_smorf_to_amp
from src.utils.log_utils import (
    get_pipeline_completion_message,
    get_task_completion_message,
)


# ============================== Custom Functions ==============================
def verify_smorf_data_availability(
    genus_name: str,
    genomes_root: str,
    logger: CustomLogger,
) -> Tuple[List[str], bool]:
    """
    Verify the availability of smorf data for a given genus.

    Parameters
    ----------
    genus_name : str
        Target genus name (e.g., "Abiotrophia").
    genomes_root : str
        Root directory containing all genus-level genome folders.
    logger : CustomLogger
        Logger instance for structured tracking.

    Returns
    -------
    (genome_dirs, success) : Tuple[List[str], bool]
        - genome_dirs : Path to the validated FASTA directory (empty if invalid)
        - success     : True if valid genome data found, False otherwise
    """
    logger.info(
        f"/ Task: Validate presence of 'smORF-derived peptides' for '{genus_name}'"
    )
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        genus_dir = os.path.join(genomes_root, genus_name)
        smorf_dir = os.path.join(genus_dir, "smorf")

        # Check if the 'smorf' subdirectory exists
        if not os.path.exists(smorf_dir):
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    f"[ smORF data summary ]\n"
                    f"▸ Genus  : '{genus_name}'\n"
                    f"▸ Reason : No genome data available for this genus"
                ),
                border="|",
                length=120,
            )
            success = False
            genome_dirs = []
        else:
            genome_dirs = sorted(
                [
                    os.path.join(smorf_dir, d)
                    for d in os.listdir(smorf_dir)
                    if os.path.isdir(os.path.join(smorf_dir, d))
                ]
            )
            n_genomes = len(genome_dirs)
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    f"[ smORF data summary ]\n"
                    f"▸ Genus  : '{genus_name}'\n"
                    f"▸ Reason : Genome data found and smORF predictions available.\n"
                    f"▸ Dirs   : {n_genomes} genome directory(ies) with smORF outputs"
                ),
                border="|",
                length=120,
            )
            success = True

        # Log completion
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        return genome_dirs, success

    except Exception:
        logger.exception("Unexpected error in 'verify_smorf_data_availability()'")
        raise


def blast_smorf_single_batch(
    genome_dirs: List[str],
    genus_name: str,
    reference_db: str,
    logger: CustomLogger,
    threads: int = 8,
    evalue: float = 5,
    word_size: int = 2,
) -> None:
    """
    Run BLASTP for all smORF-derived peptide sets under a given genus.

    Parameters
    ----------
    genome_dirs : List[str]
        List of genome-level smORF directories, each containing peptides.faa.
        Example: [..., "Abiotrophia/smorf/GCF_xxx", ...]
    genus_name : str
        Target genus name (e.g., "Abiotrophia").
    reference_db : str
        Path prefix of AMP BLAST database.
    logger : CustomLogger
        Logger instance for structured tracking.
    threads : int, optional
        Number of CPU threads (default: 8).
    evalue : float, optional
        E-value cutoff for BLASTP (default: 5).
    word_size : int, optional
        Word size for short-peptide BLASTP (default: 2).

    Returns
    -------
    None
    """
    logger.info(f"/ Task: Run BLASTP smORF to AMP for '{genus_name}'")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Summary of settings
        total = len(genome_dirs)
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"[ BLASTP Screening Summary ]\n"
                f"▸ Genus          : '{genus_name}'\n"
                f"▸ Genome count   : {total}\n"
                f"▸ Reference DB   : '{reference_db}'\n"
                f"▸ Threads        : {threads}\n"
                f"▸ E-value cutoff : {evalue}\n"
                f"▸ Word size      : {word_size}\n"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Iterate over each genome FASTA
        success_count = 0
        failed = []

        # Prepare genus-level blast root
        genus_dir = os.path.abspath(os.path.join(genome_dirs[0], "../.."))
        blast_root = os.path.join(genus_dir, "blast/smorf")
        os.makedirs(blast_root, exist_ok=True)

        # Iterate over genome-level directories
        for genome_dir in genome_dirs:
            genome_id = os.path.basename(genome_dir)
            peptides_path = os.path.join(genome_dir, "peptides.faa")

            # New: verify that peptides.faa exists
            if not os.path.exists(peptides_path):
                failed.append(genome_id)
                continue

            # Output folder: blast/<genome_id>/
            genome_blast_dir = os.path.join(blast_root, genome_id)
            os.makedirs(genome_blast_dir, exist_ok=True)
            output_csv = os.path.join(genome_blast_dir, "amp_hits.csv")

            try:

                # Run smORF for this genome
                blastp_smorf_to_amp(
                    query_fasta=peptides_path,
                    reference_db=reference_db,
                    output_csv=output_csv,
                    logger=logger,
                    threads=threads,
                    evalue=evalue,
                    word_size=word_size,
                )

                # New: verify that amp_hits.csv exists → this is the REAL success test
                if not os.path.exists(output_csv):
                    failed.append(genome_id)
                    continue

                # Otherwise count success
                success_count += 1

            except Exception:
                failed.append(genome_id)
                continue

        # Log completion summary (success / failure counts)
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"[ BLASTP Screening Completed ]\n"
                f"▸ Total genomes     : {total}\n"
                f"▸ Successfully done : {success_count}\n"
                f"▸ Failed genomes    : {len(failed)}\n"
                + (f"▸ Failed list       : {', '.join(failed)}\n" if failed else "")
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
        logger.exception("Unexpected error in 'blast_smorf_single_batch()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_smorf_to_amp_by_genus(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Run smORF-to-AMP BLAST screening pipeline for all genomes under a given genus.

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
        genus_name = kwargs.get("genomes_genus_name", "Escherichia")
        output_dir = kwargs.get(
            "genomes_output_dir",
            os.path.join(base_path, "Escherichia"),
        )
        reference_db = kwargs.get(
            "blast_reference_db",
            os.path.join(base_path, "data/processed/AMP/blast_db/amp_db"),
        )
        threads = int(kwargs.get("blast_threads", 8))
        evalue = float(kwargs.get("blast_evalue", 5))
        word_size = int(kwargs.get("blast_word_size", 2))
        genomes_root = os.path.dirname(output_dir)

        # Ensure output directory exists
        if os.path.isdir(output_dir):
            os.makedirs(name=output_dir, exist_ok=True)

        # Verify smORF data availability
        genome_dirs, has_smorf = verify_smorf_data_availability(
            genus_name=genus_name,
            genomes_root=genomes_root,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

        # Run BLAST screening (smORF → AMP)
        if has_smorf:
            blast_smorf_single_batch(
                genome_dirs=genome_dirs,
                genus_name=genus_name,
                reference_db=reference_db,
                threads=threads,
                evalue=evalue,
                word_size=word_size,
                logger=logger,
            )
            logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_smorf_to_amp_by_genus()'")
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
