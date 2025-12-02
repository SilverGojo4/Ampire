# pylint: disable=import-error, wrong-import-position, too-many-locals, consider-using-with, broad-exception-caught
"""
Download bacterial proteomes from Uniprot based on genus-level taxonomy.
"""
# ============================== Standard Library Imports ==============================
import gzip
import logging
import os
import re
import sys
import time
from collections import OrderedDict
from typing import Tuple

# ============================== Third-Party Library Imports ==============================
import pandas as pd
import requests
from ete3 import NCBITaxa

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
from src.utils.io_utils import directory_exists
from src.utils.log_utils import (
    get_pipeline_completion_message,
    get_task_completion_message,
)

# UniProt Proteome API Config
UNIPROTKB_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"


# ============================== Custom Functions ==============================
def get_genus_taxid(
    genus_name: str,
    output_path: str,
    logger: CustomLogger,
) -> Tuple[pd.DataFrame, bool]:
    """
    Retrieve the NCBI TaxID of a genus using ETE3.

    Parameters
    ----------
    genus_name : str
        Name of the genus to query (e.g., "Escherichia").
    output_path : str
        Path to save the resulting DataFrame as a CSV file.
    logger : CustomLogger
        Logger instance for progress tracking.

    Returns
    -------
    (pd.DataFrame, bool)
        A tuple (df_genus, success):
            - df_genus: DataFrame containing ['Genus', 'TaxID'] (empty if not found)
            - success : True if genus TaxID retrieved, False otherwise
    """
    logger.info("/ Task: Retrieve genus TaxID")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:

        # Initialize NCBI taxonomy database
        ncbi = NCBITaxa()

        # Retrieve genus TaxID
        name_map = ncbi.get_name_translator([genus_name])
        if not name_map or genus_name not in name_map:
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    f"Genus '{genus_name}' not found in NCBI taxonomy. "
                    f"Skipping this genus."
                ),
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            df_genus = pd.DataFrame(columns=["Genus", "TaxID"])

        else:
            genus_id = name_map[genus_name][0]

            # Create DataFrame
            df_genus = pd.DataFrame([{"Genus": genus_name, "TaxID": genus_id}])

            # Log summary
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    f"[ Taxonomy summary ]\n"
                    f"▸ Genus : '{genus_name}'\n"
                    f"▸ TaxID : {genus_id}"
                ),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

            # Ensure output directory exists
            metadata_output_dir = os.path.dirname(output_path)
            os.makedirs(metadata_output_dir, exist_ok=True)

            # Save CSV
            df_genus.to_csv(output_path, index=False)
            logger.log_with_borders(
                level=logging.INFO,
                message=f"Saved:\n'{genus_name}/metadata/genus.csv'",
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

        # Determine success flag
        success = not df_genus.empty
        return df_genus, success

    except Exception:
        logger.exception("Unexpected error in 'get_genus_taxid()'")
        raise


def split_reviewed_unreviewed(
    fasta_path: str, output_dir: str, logger: CustomLogger
) -> Tuple[bool, bool]:
    """
    Split the downloaded UniProt FASTA into reviewed (sp) and unreviewed (tr),
    optimized with streaming writes to avoid storing large buffers in memory.

    Parameters
    ----------
    fasta_path : str
        Path to the downloaded UniProt FASTA file (typically 'fasta/all_raw.fasta').
    output_dir : str
        Directory where the split FASTA files will be written.
    logger : CustomLogger
        Logger instance for progress tracking.

    Returns
    -------
    (bool, bool)
        has_reviewed, has_unreviewed
        - True if this category contains >=1 entry
        - False if empty (no such sequences)
    """
    try:
        reviewed_path = os.path.join(output_dir, "reviewed.fasta")
        unreviewed_path = os.path.join(output_dir, "unreviewed.fasta")

        reviewed_count = 0
        unreviewed_count = 0

        current_buffer = []
        current_is_reviewed = None

        # Choose correct opener: gzip or normal text
        def open_fasta(path: str):
            if path.endswith(".gz"):
                return gzip.open(path, "rt", encoding="utf-8")
            return open(path, "r", encoding="utf-8")

        # Use context managers for file writers (pylint R1732)
        with open(reviewed_path, "w", encoding="utf-8") as rf, open(
            unreviewed_path, "w", encoding="utf-8"
        ) as uf:

            def flush_entry():
                nonlocal reviewed_count, unreviewed_count
                if not current_buffer:
                    return

                entry_text = "\n".join(current_buffer) + "\n"

                if current_is_reviewed:
                    rf.write(entry_text)
                    reviewed_count += 1
                else:
                    uf.write(entry_text)
                    unreviewed_count += 1

            # Parse FASTA entries (streaming, memory-efficient)
            with open_fasta(fasta_path) as f:
                for line in f:
                    line = line.rstrip("\n")

                    if line.startswith(">"):
                        # flush previous entry
                        flush_entry()

                        # start new entry
                        current_buffer = [line]
                        current_is_reviewed = line.startswith(">sp|")

                    else:
                        current_buffer.append(line)

                # flush last entry at EOF
                flush_entry()

        # Create subdirectories for reviewed / unreviewed
        has_reviewed = reviewed_count > 0
        has_unreviewed = unreviewed_count > 0

        if has_reviewed:
            os.makedirs(os.path.join(output_dir, "reviewed"), exist_ok=True)
        if has_unreviewed:
            os.makedirs(os.path.join(output_dir, "unreviewed"), exist_ok=True)

        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Split Summary ]\n"
                f"▸ Reviewed sequences   : {reviewed_count}\n"
                f"▸ Unreviewed sequences : {unreviewed_count}"
            ),
            border="|",
            length=120,
        )

        return has_reviewed, has_unreviewed

    except Exception:
        logger.exception("Unexpected error in 'split_reviewed_unreviewed()'")
        raise


def split_fasta_by_taxid(
    fasta_path: str,
    output_dir: str,
    logger: CustomLogger,
) -> None:
    """
    Split a UniProt FASTA file (reviewed.fasta or unreviewed.fasta)
    into multiple FASTA files grouped by NCBI TaxID (OX=), using streaming
    writes to avoid large buffer accumulation in memory.

    Parameters
    ----------
    fasta_path : str
        Path to the input UniProt FASTA file.
    output_dir : str
        Directory where species-level FASTA files will be written.
    logger : CustomLogger
        Logger instance for progress tracking.

    Returns
    -------
    None
    """
    try:
        # LRU file handle cache: taxid → open file handle
        file_handles = OrderedDict()
        max_open_files = 64  # prevent "Too many open files"
        species_count = 0

        # current entry state
        current_taxid = None
        current_buffer = []

        def get_handle(taxid: str):
            """
            Retrieve file handle using LRU strategy.
            If too many files are open, close the least-recently-used one.
            """
            nonlocal species_count

            # Case 1 — already open → mark as recently used
            if taxid in file_handles:
                fh = file_handles.pop(taxid)
                file_handles[taxid] = fh
                return fh

            # Case 2 — need to open a new file
            if len(file_handles) >= max_open_files:
                _, old_fh = file_handles.popitem(last=False)
                old_fh.close()

            out_path = os.path.join(output_dir, f"{taxid}.fasta")

            # pylint fix: pre-create file using with-statement
            with open(out_path, "a", encoding="utf-8"):
                pass

            # reopen for persistent handle (cannot use with here)
            fh = open(out_path, "a", encoding="utf-8")
            file_handles[taxid] = fh
            species_count += 1
            return fh

        # Parse FASTA entries
        with open(fasta_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")

                # New entry
                if line.startswith(">"):

                    # flush previous entry
                    if current_taxid is not None and current_buffer:
                        fh = get_handle(current_taxid)
                        fh.write("\n".join(current_buffer))
                        fh.write("\n")

                    # start new entry buffer
                    current_buffer = [line]

                    # extract TaxID (OX=XXXX)
                    taxid_match = re.search(r"OX=(\d+)", line)
                    current_taxid = taxid_match.group(1) if taxid_match else None

                else:
                    current_buffer.append(line)

            # flush last entry
            if current_taxid is not None and current_buffer:
                fh = get_handle(current_taxid)
                fh.write("\n".join(current_buffer))
                fh.write("\n")

        # close all open file handles
        for fh in file_handles.values():
            fh.close()

        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Species Split by TaxID ]\n"
                f"▸ Input FASTA      : '{fasta_path}'\n"
                f"▸ Output directory : '{output_dir}'\n"
                f"▸ TaxIDs (species) : {species_count}\n"
                f"▸ Max open files   : {max_open_files}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error in 'split_fasta_by_taxid()'")
        raise


def cleanup_intermediate_fastas(
    fasta_dir: str,
    keep_raw: bool = False,
    keep_split: bool = False,
) -> None:
    """
    Clean up intermediate FASTA files created during UniProt proteome downloading.

    Parameters
    ----------
    fasta_dir : str
        Directory where FASTA files are stored (typically <genus>/fasta/)
    keep_raw : bool
        If True, keep all_raw.fasta
    keep_split : bool
        If True, keep reviewed.fasta and unreviewed.fasta

    Returns
    -------
    None
    """

    # RAW group
    raw_files = ["all_raw.fasta", "all_raw.fasta.gz"]

    # SPLIT group
    split_files = ["reviewed.fasta", "unreviewed.fasta"]

    # Decide which to delete
    files_to_delete = []

    if not keep_raw:
        files_to_delete.extend(raw_files)

    if not keep_split:
        files_to_delete.extend(split_files)

    # Perform deletion
    for fname in files_to_delete:
        fpath = os.path.join(fasta_dir, fname)
        if os.path.exists(fpath):
            try:
                os.remove(fpath)
            except Exception:
                pass


def uniprot_proteome_download(
    genus_taxid: int,
    output_dir: str,
    logger: CustomLogger,
    include_isoform: bool = False,
    chunk_size: int = 262144,
) -> None:
    """
    Download ALL UniProtKB protein sequences (reviewed + unreviewed)
    for a genus-level TaxID using a SINGLE UniProt stream query.

    Parameters
    ----------
    genus_taxid : int
        NCBI TaxID of the genus.
    output_dir : str
        Directory where proteome files and metadata will be stored.
    logger : CustomLogger
        Project-specific logger for structured logging.

    Returns
    -------
    None
    """
    logger.info("/ Task: Run 'uniprot-proteome-download' for bacterial proteomes")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Ensure output directory exists
        fasta_output_dir = os.path.join(output_dir, "fasta")
        metadata_output_dir = os.path.join(output_dir, "metadata")
        for sub_output_dir in [fasta_output_dir, metadata_output_dir]:
            os.makedirs(sub_output_dir, exist_ok=True)

        # Output FASTA path (gzipped)
        fasta_path = os.path.join(fasta_output_dir, "all_raw.fasta.gz")

        # Log parameters and command
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ UniProt Proteome Query Parameters ]\n"
                f"▸ TaxID : {genus_taxid}\n"
                f"▸ Include isoforms : '{include_isoform}'\n"
                f"▸ Compressed (gzip) : true"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Prepare API query
        params = {
            "query": f"taxonomy_id:{genus_taxid}",
            "format": "fasta",
            "compressed": "true",
        }
        if include_isoform:
            params["includeIsoform"] = "true"

        headers = {
            "User-Agent": "Ampire-Uniprot-Downloader (contact: your_email@example.com)",
            "Accept-Encoding": "gzip",
        }

        # Download stream
        stream_received = False
        with requests.get(
            UNIPROTKB_STREAM_URL,
            params=params,
            headers=headers,
            stream=True,
            timeout=300,
        ) as r:

            # HTTP status check
            try:
                r.raise_for_status()
            except Exception as exc:
                snippet = r.content[:1000]
                raise RuntimeError(
                    f"HTTP Error (gzipped response):\n{snippet!r}"
                ) from exc

            with open(fasta_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        stream_received = True
                        f.write(chunk)

        if not stream_received:
            logger.log_with_borders(
                level=logging.INFO,
                message=(
                    f"No UniProt protein sequences found for TaxID={genus_taxid}. "
                    "This genus may not have entries in Uniprot."
                ),
                border="|",
                length=120,
            )

        else:
            # Check file existence
            if not os.path.exists(fasta_path):
                raise RuntimeError(
                    f"UniProt returned data but FASTA file was NOT written: '{fasta_path}'"
                )

            # Check file size
            file_size = os.path.getsize(fasta_path)
            if file_size == 0:
                raise RuntimeError(
                    f"Gzipped FASTA file exists but is EMPTY:  '{fasta_path}'"
                )

            # Split into reviewed / unreviewed
            has_reviewed, has_unreviewed = split_reviewed_unreviewed(
                fasta_path=fasta_path,
                output_dir=fasta_output_dir,
                logger=logger,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            if has_reviewed:
                reviewed_path = os.path.join(fasta_output_dir, "reviewed.fasta")
                reviewed_dir = os.path.join(fasta_output_dir, "reviewed")
                split_fasta_by_taxid(
                    fasta_path=reviewed_path, output_dir=reviewed_dir, logger=logger
                )
            if has_unreviewed:
                unreviewed_path = os.path.join(fasta_output_dir, "unreviewed.fasta")
                unreviewed_dir = os.path.join(fasta_output_dir, "unreviewed")
                split_fasta_by_taxid(
                    fasta_path=unreviewed_path, output_dir=unreviewed_dir, logger=logger
                )

        # Log completion
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error during 'uniprot_proteome_download()'.")
        raise


# ============================== Pipeline Entry Point ==============================
def run_download_proteomes_by_genus(
    base_path: str, logger: CustomLogger, **kwargs
) -> None:
    """
    Download all bacterial proteomes from Uniprot for a specified genus.

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
        keep_raw = kwargs.get("keep_raw", False)
        keep_split = kwargs.get("keep_split", False)

        # Ensure output directory exists
        if not directory_exists(dir_path=output_dir):
            os.makedirs(name=output_dir, exist_ok=True)

        # Retrieve species under the genus
        genus_csv = os.path.join(output_dir, "metadata/genus.csv")
        df_genus, success = get_genus_taxid(genus_name, genus_csv, logger)
        logger.add_spacer(level=logging.INFO, lines=1)

        # Download proteomes
        if success:
            genus_taxid = int(df_genus["TaxID"].iloc[0])
            uniprot_proteome_download(genus_taxid, output_dir, logger)

            # cleanup intermediate files
            fasta_dir = os.path.join(output_dir, "fasta")
            cleanup_intermediate_fastas(
                fasta_dir=fasta_dir,
                keep_raw=keep_raw,
                keep_split=keep_split,
            )
            logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_download_proteomes_by_genus()'")
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
