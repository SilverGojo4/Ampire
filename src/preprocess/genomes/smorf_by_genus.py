# pylint: disable=import-error, wrong-import-position, line-too-long, too-many-arguments, too-many-positional-arguments, too-many-locals, too-many-statements
"""
Predict small open reading frames (smORFs) from bacterial genomes on a genus-wide scale.
"""
# ============================== Standard Library Imports ==============================
import csv
import glob
import logging
import os
import shutil
import subprocess
import sys
import time
from typing import Dict, List, Tuple

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
from src.utils.log_utils import (
    get_pipeline_completion_message,
    get_task_completion_message,
)


# ============================== Custom Functions ==============================
def verify_genome_data_availability(
    genus_name: str,
    genomes_root: str,
    logger: CustomLogger,
) -> Tuple[str, bool]:
    """
    Verify the availability of genome data for a given genus.

    This function checks whether:
        1. The genus directory exists under `genomes_root`.
        2. The "fasta" subdirectory exists (implying that NCBI returned results).
        3. The "fasta" subdirectory contains genome FASTA files (indicating available assemblies).

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
    (fasta_dir, success) : Tuple[str, bool]
        - fasta_dir : Path to the validated FASTA directory (empty if invalid)
        - success   : True if valid genome data found, False otherwise
    """
    logger.info(f"/ Task: Validate existence of 'genome data' for '{genus_name}'")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        genus_dir = os.path.join(genomes_root, genus_name)
        fasta_dir = os.path.join(genus_dir, "fasta")

        # Check if the 'fasta' subdirectory exists
        if not os.path.exists(fasta_dir):
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    f"[ Genome data summary ]\n"
                    f"▸ Genus  : '{genus_name}'\n"
                    f"▸ Reason : This genus could not be found in NCBI taxonomy"
                ),
                border="|",
                length=120,
            )
            success = False
            fasta_dir = ""
        else:
            fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith((".fna"))]
            if not fasta_files:
                logger.log_with_borders(
                    level=logging.INFO,
                    message=(
                        f"[ Genome data summary ]\n"
                        f"▸ Genus  : '{genus_name}'\n"
                        f"▸ Reason : Species were found in NCBI, but no genome assemblies were available or retrieved."
                    ),
                    border="|",
                    length=120,
                )
                success = False
                fasta_dir = ""
            else:
                n_files = len(fasta_files)
                logger.log_with_borders(
                    level=logging.INFO,
                    message=(
                        f"[ Genome data summary ]\n"
                        f"▸ Genus  : '{genus_name}'\n"
                        f"▸ Reason : Genome data successfully retrieved and validated.\n"
                        f"▸ Files  : {n_files} genome FASTA file(s)"
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

        return fasta_dir, success

    except Exception:
        logger.exception("Unexpected error in 'verify_genome_data_availability()'")
        raise


def smorf_single_genome(
    fasta_path: str,
    output_dir: str,
    logger: CustomLogger,
    dsn1_cutoff: float = 0.9999,
    dsn2_cutoff: float = 0.9999,
    phmm_cutoff: float = 1e-6,
    overwrite: bool = False,
) -> bool:
    """
    Run smORF on a single genome FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to the input genome FASTA file (.fna, .fa, .fasta).
    output_dir : str
        Directory where smORF output will be stored.
    logger : CustomLogger
        Logger instance for structured tracking.
    dsn1_cutoff, dsn2_cutoff : float, optional
        Probability cutoffs for DSN1/DSN2 models.
    phmm_cutoff : float, optional
        Profile HMM significance cutoff.
    overwrite : bool, optional
        Whether to overwrite existing outputs.

    Returns
    -------
    bool
        True if execution succeeded, False otherwise.
    """
    try:

        # Build SmORF command
        cmd = [
            "smorf",
            "single",
            fasta_path,
            "-o",
            output_dir,
            "-idsn1",
            str(dsn1_cutoff),
            "-idsn2",
            str(dsn2_cutoff),
            "-iphmm",
            str(phmm_cutoff),
        ]
        if overwrite:
            cmd.append("--force")

        # Run subprocess
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        # Explicit success/failure handling
        if result.returncode == 0:
            return True

        return False

    except Exception:
        logger.exception(
            f"Unexpected error during 'smorf_single_genome()' for '{os.path.splitext(os.path.basename(fasta_path))[0]}'"
        )
        raise


def cleanup_smorf_tmp(output_dir: str) -> None:
    """
    Clean up unnecessary intermediate files in the smORF 'tmp/' directory.

    Parameters
    ----------
    output_dir : str
        Path to the genome-level smORF output directory.
    """
    tmp_dir = os.path.join(output_dir, "tmp")
    remove_targets = [
        "prodigal.faa",
        "prodigal.ffn",
        "prodigal.gff",
    ]

    for fname in remove_targets:
        fpath = os.path.join(tmp_dir, fname)
        if os.path.exists(fpath):
            os.remove(fpath)


def rename_smorf_outputs(
    output_dir: str,
) -> None:
    """
    Rename core smORF output files to standardized names.

    Parameters
    ----------
    output_dir : str
        Path to the genome-level smORF output directory.

    Returns
    -------
    None
    """
    # Define renaming rules
    rename_map = {
        ".faa": "peptides.faa",
        ".ffn": "nucleotides.ffn",
        ".gff": "annotations.gff",
        ".tsv": "predictions.tsv",
    }

    # Apply renaming
    for ext, new_name in rename_map.items():
        matches = glob.glob(os.path.join(output_dir, f"*{ext}"))
        if not matches:
            continue
        old_path = matches[0]
        new_path = os.path.join(output_dir, new_name)

        # Avoid overwriting if target already exists
        if os.path.exists(new_path):
            continue

        shutil.move(old_path, new_path)


def read_tsv(
    path: str,
    delimiter: str = "\t",
    comment: str = "#",
) -> List[Dict[str, str]]:
    """
    Read a tab-separated file into a list of dictionaries.

    Parameters
    ----------
    path : str
        Path to the input TSV file.
    delimiter : str, optional
        Field delimiter used in the file (default: '\\t').
    comment : str, optional
        Comment character; lines beginning with this character are ignored (default: '#').

    Returns
    -------
    List[Dict[str, str]]
        List of records, each represented as a dictionary where keys are column indices (as strings)
        and values are corresponding cell contents.
    """
    records = []

    # Open and parse file line-by-line
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(comment) or not line.strip():
                continue
            parts = line.strip().split(delimiter)
            # Store as dictionary using column indices as keys
            records.append({str(i): val for i, val in enumerate(parts)})

    return records


def save_smorf_table(
    records: List[Dict[str, str]],
    out_path: str,
) -> None:
    """
    Write a list of smORF prediction records to a TSV file.

    Parameters
    ----------
    records : List[Dict[str, str]]
        List of smORF prediction records to write.
    out_path : str
        Path to the output TSV file.

    Returns
    -------
    None
    """
    if not records:
        return

    # Ensure output directory exists
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # Write header and rows
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def assemble_smorf_predictions_no_pandas(
    genome_dir: str,
) -> List[Dict[str, str]]:
    """
    Assemble the final smORF predictions table by merging:
    - tmp/model_predictions.tsv
    - tmp/prodigal.small.gff
    - tmp/hmmsearch.tbl

    Parameters
    ----------
    genome_dir : str
        Path to the genome-level smORF output directory (e.g., "Abiotrophia/smorf/GCF_xxx").

    Returns
    -------
    List[Dict[str, str]]
        Combined list of smORF prediction records containing metadata, coordinates,
        model probabilities, and HMM annotations.
    """
    tmp_dir = os.path.join(genome_dir, "tmp")

    # Load model predictions
    model_path = os.path.join(tmp_dir, "model_predictions.tsv")
    with open(model_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        model_records = [dict(r) for r in reader]

    for r in model_records:
        if "orf_seq" in r:
            r["orf"] = r.pop("orf_seq")  # Standardize field name

    # Parse prodigal.small.gff
    gff_path = os.path.join(tmp_dir, "prodigal.small.gff")
    gff_records = []

    with open(gff_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue

            seqid = ""
            attrs = cols[8]
            if "ID=" in attrs:
                seqid = attrs.split("ID=")[1].split(";")[0]

            gff_records.append(
                {
                    "seqid": seqid,
                    "contig": cols[0],
                    "start": cols[3],
                    "end": cols[4],
                    "orient": cols[6],
                }
            )

    # Parse hmmsearch.tbl
    hmm_path = os.path.join(tmp_dir, "hmmsearch.tbl")
    hmm_records = {}

    with open(hmm_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split()
            seqid, smorfam, evalue = cols[0], cols[2], cols[4]

            # Keep the lowest e-value per sequence ID
            if seqid not in hmm_records or float(evalue) < float(
                hmm_records[seqid]["hmm_smorfam_evalue"]
            ):
                hmm_records[seqid] = {
                    "smorfam": smorfam,
                    "hmm_smorfam_evalue": evalue,
                }

    # Merge all sources
    merged = []
    gff_dict = {r["seqid"]: r for r in gff_records}

    for rec in model_records:
        seqid = rec.get("seqid", "")
        merged_rec = {
            "seqid": seqid,
            "contig": gff_dict.get(seqid, {}).get("contig", ""),
            "start": gff_dict.get(seqid, {}).get("start", ""),
            "end": gff_dict.get(seqid, {}).get("end", ""),
            "orient": gff_dict.get(seqid, {}).get("orient", ""),
            "smorfam": hmm_records.get(seqid, {}).get("smorfam", ""),
            "hmm_smorfam_evalue": hmm_records.get(seqid, {}).get(
                "hmm_smorfam_evalue", ""
            ),
            "dsn1_prob_smorf": rec.get("dsn1_prob_smorf", ""),
            "dsn2_prob_smorf": rec.get("dsn2_prob_smorf", ""),
            "5p_seq": rec.get("5p_seq", ""),
            "orf": rec.get("orf", ""),
            "3p_seq": rec.get("3p_seq", ""),
        }
        merged.append(merged_rec)

    return merged


def smorf_single_batch(
    fasta_dir: str,
    genus_name: str,
    genomes_root: str,
    logger: CustomLogger,
    dsn1_cutoff: float = 0.9999,
    dsn2_cutoff: float = 0.9999,
    phmm_cutoff: float = 1e-6,
    overwrite: bool = False,
) -> None:
    """
    Iterate over all genome FASTA files under `fasta_dir`
    and list their names (preparation for smORF batch execution).

    Parameters
    ----------
    fasta_dir : str
        Directory containing genome FASTA files (.fna, .fa, .fasta).
    genus_name : str
        Target genus name (e.g., "Abiotrophia").
    genomes_root : str
        Root directory containing all genus-level genome folders.
    logger : CustomLogger
        Logger instance for structured tracking.
    dsn1_cutoff : float, optional
        DSN1 individual model probability cutoff (default: 0.9999).
    dsn2_cutoff : float, optional
        DSN2 individual model probability cutoff (default: 0.9999).
    phmm_cutoff : float, optional
        Profile HMM significance cutoff (default: 1e-6).
    overwrite : bool, optional
        Whether to overwrite existing output directory (default: False).

    Returns
    -------
    None
    """
    logger.info(f"/ Task: Run smORF for '{genus_name}'")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Ensure output directories exist
        genus_dir = os.path.join(genomes_root, genus_name)
        smorf_dir = os.path.join(genus_dir, "smorf")
        smorf_all_dir = os.path.join(genus_dir, "smorf_all")
        for dir_path in [smorf_dir, smorf_all_dir]:
            os.makedirs(name=dir_path, exist_ok=True)

        # Collect all FASTA files
        fasta_files = sorted([f for f in os.listdir(fasta_dir) if f.endswith(".fna")])
        total = len(fasta_files)

        # Genus-level summary
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"[ smORF Execution Summary ]\n"
                f"▸ Genus        : '{genus_name}'\n"
                f"▸ Genome count : {total}\n"
                f"▸ DSN1 cutoff  : {dsn1_cutoff}\n"
                f"▸ DSN2 cutoff  : {dsn2_cutoff}\n"
                f"▸ pHMM cutoff  : {phmm_cutoff}\n"
                f"▸ Overwrite    : {overwrite}\n"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Iterate over each genome FASTA
        success_count = 0
        failed = []

        # Iterate over each genome FASTA
        for fasta_file in fasta_files:
            genome_id = os.path.splitext(fasta_file)[0]
            fasta_path = os.path.join(fasta_dir, fasta_file)
            smorf_output_dir = os.path.join(smorf_dir, genome_id)
            smorf_all_output_dir = os.path.join(smorf_all_dir, genome_id)

            # Run smORF for this genome
            _ = smorf_single_genome(
                fasta_path=fasta_path,
                output_dir=smorf_output_dir,
                logger=logger,
                dsn1_cutoff=dsn1_cutoff,
                dsn2_cutoff=dsn2_cutoff,
                phmm_cutoff=phmm_cutoff,
                overwrite=overwrite,
            )

            # New: verify that model_predictions.tsv exists → this is the REAL success test
            model_file = os.path.join(smorf_output_dir, "tmp", "model_predictions.tsv")
            if not os.path.exists(model_file):
                failed.append(genome_id)
                continue

            # Otherwise count success
            success_count += 1

            # Clean up unnecessary intermediate files
            cleanup_smorf_tmp(output_dir=smorf_output_dir)

            # Assemble final prediction table
            records = assemble_smorf_predictions_no_pandas(genome_dir=smorf_output_dir)
            output_path = os.path.join(smorf_all_output_dir, "predictions.tsv")
            save_smorf_table(records=records, out_path=output_path)

            # Copy key small-ORF files from tmp to smorf_all/
            for fname in [
                "prodigal.small.faa",
                "prodigal.small.ffn",
                "prodigal.small.gff",
            ]:
                src = os.path.join(smorf_output_dir, "tmp", fname)
                dst = os.path.join(smorf_all_output_dir, fname)
                if os.path.exists(src):
                    shutil.copy(src, dst)

            # Remove the entire tmp/ directory for cleanliness
            shutil.rmtree(os.path.join(smorf_output_dir, "tmp"), ignore_errors=True)

            # Rename outputs in both directories
            rename_smorf_outputs(output_dir=smorf_output_dir)
            rename_smorf_outputs(output_dir=smorf_all_output_dir)

            # Convert predictions.tsv → predictions.csv in both directories
            for dir_path in [smorf_output_dir, smorf_all_output_dir]:
                tsv_path = os.path.join(dir_path, "predictions.tsv")
                csv_path = os.path.join(dir_path, "predictions.csv")

                if os.path.exists(tsv_path):
                    # Read TSV manually
                    with open(tsv_path, "r", newline="", encoding="utf-8") as f_in:
                        reader = csv.reader(f_in, delimiter="\t")
                        rows = list(reader)

                    # Write as CSV (comma-separated)
                    with open(csv_path, "w", newline="", encoding="utf-8") as f_out:
                        writer = csv.writer(f_out, delimiter=",")
                        writer.writerows(rows)

                    # Optionally remove old .tsv to avoid redundancy
                    os.remove(tsv_path)

        # Log completion summary (success / failure counts)
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"[ smORF Completed ]\n"
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
        logger.exception("Unexpected error in 'smorf_single_batch()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_smorf_by_genus(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Run SmORFinder prediction pipeline for all genomes under a specified genus.

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
        genomes_root = os.path.dirname(output_dir)

        # Ensure output directory exists
        if os.path.isdir(output_dir):
            os.makedirs(name=output_dir, exist_ok=True)

        # Verify genome data availability
        fasta_dir, has_genomes = verify_genome_data_availability(
            genus_name=genus_name,
            genomes_root=genomes_root,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

        # Run smORF
        if has_genomes:
            smorf_single_batch(
                fasta_dir=fasta_dir,
                genus_name=genus_name,
                genomes_root=genomes_root,
                logger=logger,
            )
            logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_smorf_by_genus()'")
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
