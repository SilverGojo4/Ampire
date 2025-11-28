# pylint: disable=import-error, wrong-import-position, too-many-arguments, too-many-positional-arguments
"""
Blastp Runner Module
"""
# ============================== Standard Library Imports ==============================
import os
import sys
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


# ============================== Custom Functions ==============================
def blastp_smorf_to_amp(
    query_fasta: str,
    reference_db: str,
    output_csv: str,
    logger: CustomLogger,
    threads: int = 8,
    evalue: float = 5,
    word_size: int = 2,
) -> None:
    """
    Identify AMP homologs from smORF-derived peptides using a short-peptide-
    optimized BLASTP workflow.

    Parameters
    ----------
    query_fasta : str
        Path to smORF peptides FASTA (query sequences).
    reference_db : str
        Path prefix of the AMP BLAST database (without file extension).
    output_csv : str
        Path to save the parsed BLASTP results as CSV.
    logger : CustomLogger
        Logger instance for structured tracking.
    threads : int
        Number of CPU threads.
    evalue : float
        BLAST E-value cutoff (default=5), recommended for short-peptide
        sensitivity.
    word_size : int
        Word size for BLASTP seed matching (default=2, required for short peptides).

    Returns
    -------
    None
    """
    try:
        # Prepare tmp TSV path
        output_dir = os.path.dirname(output_csv)
        os.makedirs(output_dir, exist_ok=True)
        tmp_tsv = output_csv.replace(".csv", ".tsv")

        # Run BLASTP
        cmd = [
            "blastp",
            "-query",
            query_fasta,
            "-db",
            reference_db,
            "-out",
            tmp_tsv,
            "-num_threads",
            str(threads),
            "-evalue",
            str(evalue),
            "-word_size",
            str(word_size),
            "-matrix",
            "PAM30",
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen "
            "qstart qend sstart send evalue bitscore qlen slen",
        ]
        subprocess.run(cmd, check=True)

        # Read TSV without header
        column_names = [
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qlen",
            "slen",
        ]
        df = pd.read_csv(tmp_tsv, sep="\t", header=None, names=column_names)

        # Save to final CSV
        df.to_csv(output_csv, index=False)

        # Remove tmp TSV
        os.remove(tmp_tsv)

    except Exception:
        logger.exception(
            f"Unexpected error during 'blastp_smorf_to_amp()' for '{query_fasta}'"
        )
        raise
