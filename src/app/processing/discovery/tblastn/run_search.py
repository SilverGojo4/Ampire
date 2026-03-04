# pylint: disable=too-many-arguments, too-many-locals
"""
TBLASTN discovery search utilities.

This module implements the production-stage logic for running a tblastn
search against a pre-built BLAST nucleotide database, and writing the
raw tabular output to disk.
"""

from __future__ import annotations

# Standard Library Imports
import shlex
import subprocess
from pathlib import Path
from typing import List, Sequence

# Default outfmt 6 fields emitted by tblastn discovery searches
DEFAULT_OUTFMT_FIELDS: List[str] = [
    "qseqid",
    "qlen",
    "sseqid",
    "slen",
    "sframe",
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
    "qcovhsp",
    "qseq",
    "sseq",
]


# Public API
def run_tblastn_search(
    *,
    query_fasta: Path,
    db_path: Path,
    out_tsv: Path,
    task: str = "tblastn",
    evalue: float = 100000,
    word_size: int = 2,
    seg: str = "no",
    soft_masking: bool = False,
    comp_based_stats: int = 0,
    use_sw_tback: bool = True,
    max_target_seqs: int = 200,
    max_hsps: int = 50,
    qcov_hsp_perc: float = 0,
    num_threads: int = 32,
    outfmt_fields: Sequence[str] | None = None,
    overwrite: bool = False,
) -> Path:
    """
    Execute a tblastn discovery search and write raw tabular output.

    Parameters
    ----------
    query_fasta : pathlib.Path
        Protein query FASTA file used as tblastn input.
    db_path : pathlib.Path
        BLAST database prefix path (e.g., /path/to/db_prefix).
    out_tsv : pathlib.Path
        Output TSV file written using BLAST outfmt 6.
    task : str, optional
        BLAST task mode (default: "tblastn").
    evalue : float, optional
        E-value threshold applied during search.
    word_size : int, optional
        Word size parameter controlling seed matching.
    seg : str, optional
        SEG filtering mode (e.g., "no").
    soft_masking : bool, optional
        Whether soft masking is enabled.
    comp_based_stats : int, optional
        Composition-based statistics mode.
    use_sw_tback : bool, optional
        Enable Smith–Waterman traceback refinement.
    max_target_seqs : int, optional
        Maximum target sequences reported per query.
    max_hsps : int, optional
        Maximum HSPs reported per subject sequence.
    num_threads : int, optional
        Number of CPU threads used by tblastn.
    outfmt_fields : Sequence[str], optional
        Fields included in BLAST outfmt 6 output.
    overwrite : bool, optional
        Whether to overwrite an existing output file.

    Returns
    -------
    pathlib.Path
        Path to the generated raw alignment TSV file.

    Raises
    ------
    FileNotFoundError
        If query FASTA or BLAST database prefix is invalid.
    FileExistsError
        If output exists and overwrite is False.
    ValueError
        If provided parameters are invalid.
    RuntimeError
        If tblastn execution fails.
    """
    if outfmt_fields is None:
        outfmt_fields = DEFAULT_OUTFMT_FIELDS

    query_fasta = Path(query_fasta)
    db_path = Path(db_path)
    out_tsv = Path(out_tsv)

    _validate_inputs(query_fasta=query_fasta, db_path=db_path)
    _validate_params(
        task=task,
        evalue=evalue,
        word_size=word_size,
        comp_based_stats=comp_based_stats,
        max_target_seqs=max_target_seqs,
        max_hsps=max_hsps,
        num_threads=num_threads,
        outfmt_fields=outfmt_fields,
    )

    # Prepare output location
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    if out_tsv.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists and overwrite is disabled: '{out_tsv}'"
        )

    if out_tsv.exists() and overwrite:
        out_tsv.unlink(missing_ok=True)

    # Build command
    cmd = _build_tblastn_cmd(
        query_fasta=query_fasta,
        db_path=db_path,
        out_tsv=out_tsv,
        task=task,
        evalue=evalue,
        word_size=word_size,
        seg=seg,
        soft_masking=soft_masking,
        comp_based_stats=comp_based_stats,
        use_sw_tback=use_sw_tback,
        max_target_seqs=max_target_seqs,
        max_hsps=max_hsps,
        qcov_hsp_perc=qcov_hsp_perc,
        num_threads=num_threads,
        outfmt_fields=outfmt_fields,
    )

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )

    stderr = (result.stderr or "").strip()
    if result.returncode != 0:
        raise RuntimeError(
            "tblastn failed.\n"
            f"Command: {shlex.join(cmd)}\n"
            f"Return code: {result.returncode}\n"
            f"STDERR: {stderr}"
        )

    # Empty output is allowed (no hits scenario)
    return out_tsv


# Internal Helpers
def _validate_inputs(*, query_fasta: Path, db_path: Path) -> None:
    """
    Validate existence of query FASTA and BLAST database prefix.

    Parameters
    ----------
    query_fasta : pathlib.Path
        Protein query FASTA file used as tblastn input.
    db_path : pathlib.Path
        BLAST database prefix path. The function checks for
        existence of associated BLAST database index files.

    Raises
    ------
    FileNotFoundError
        If the query FASTA file does not exist, the database
        directory is missing, or no valid BLAST database files
        matching the prefix are detected.
    """

    # Validate query FASTA
    if not query_fasta.is_file():
        raise FileNotFoundError(f"Query FASTA not found: '{query_fasta}'")

    # Validate database directory
    db_dir = db_path.parent
    db_name = db_path.name

    if not db_dir.exists():
        raise FileNotFoundError(f"DB directory not found: '{db_dir}'")

    # Check BLAST alias files (.nal / .pal)
    nal_file = db_dir / f"{db_name}.nal"
    pal_file = db_dir / f"{db_name}.pal"

    if nal_file.exists() or pal_file.exists():
        return

    # Check volume-split BLAST DB files
    # (e.g., db.00.nin, db.01.nsq, ...)
    required_ext = [".nin", ".nsq", ".nhr"]

    for ext in required_ext:
        matches = list(db_dir.glob(f"{db_name}.*{ext}"))
        if not matches:
            raise FileNotFoundError(
                "Missing BLAST DB index files.\n"
                f"Expected files matching pattern: '{db_name}.*{ext}'\n"
                f"DB prefix: '{db_path}'"
            )


def _validate_params(
    *,
    task: str,
    evalue: float,
    word_size: int,
    comp_based_stats: int,
    max_target_seqs: int,
    max_hsps: int,
    num_threads: int,
    outfmt_fields: Sequence[str],
) -> None:
    """
    Validate tblastn execution parameters.

    Parameters
    ----------
    task : str
        BLAST task name.
    evalue : float
        E-value threshold (> 0).
    word_size : int
        Seed word size (> 0).
    comp_based_stats : int
        Composition-based statistics mode (>= 0).
    max_target_seqs : int
        Maximum reported target sequences (> 0).
    max_hsps : int
        Maximum HSPs per subject (> 0).
    num_threads : int
        Number of execution threads (> 0).
    outfmt_fields : Sequence[str]
        Fields defining BLAST outfmt 6 schema.

    Raises
    ------
    ValueError
        If any parameter violates required constraints.
    """
    if not task:
        raise ValueError("task must be a non-empty string.")
    if evalue <= 0:
        raise ValueError("evalue must be > 0.")
    if word_size <= 0:
        raise ValueError("word_size must be > 0.")
    if comp_based_stats < 0:
        raise ValueError("comp_based_stats must be >= 0.")
    if max_target_seqs <= 0:
        raise ValueError("max_target_seqs must be > 0.")
    if max_hsps <= 0:
        raise ValueError("max_hsps must be > 0.")
    if num_threads <= 0:
        raise ValueError("num_threads must be > 0.")
    if not outfmt_fields:
        raise ValueError("outfmt_fields must not be empty.")


def _format_outfmt(fields: Sequence[str]) -> str:
    """
    Format BLAST outfmt 6 specification string.

    Parameters
    ----------
    fields : Sequence[str]
        Ordered list of BLAST output fields.

    Returns
    -------
    str
        Formatted outfmt string compatible with tblastn.
    """
    return "6 " + " ".join(fields)


def _build_tblastn_cmd(
    *,
    query_fasta: Path,
    db_path: Path,
    out_tsv: Path,
    task: str,
    evalue: float,
    word_size: int,
    seg: str,
    soft_masking: bool,
    comp_based_stats: int,
    use_sw_tback: bool,
    max_target_seqs: int,
    max_hsps: int,
    qcov_hsp_perc: float,
    num_threads: int,
    outfmt_fields: Sequence[str],
) -> List[str]:
    """
    Construct tblastn command-line arguments.

    Parameters
    ----------
    query_fasta : pathlib.Path
        Protein query FASTA file.
    db_path : pathlib.Path
        BLAST database prefix.
    out_tsv : pathlib.Path
        Output TSV file path.
    task : str
        BLAST task mode.
    evalue : float
        E-value threshold.
    word_size : int
        Word size parameter.
    seg : str
        SEG filtering configuration.
    soft_masking : bool
        Whether soft masking is enabled.
    comp_based_stats : int
        Composition-based statistics mode.
    use_sw_tback : bool
        Enable Smith–Waterman traceback.
    max_target_seqs : int
        Maximum target sequences reported.
    max_hsps : int
        Maximum HSPs per subject.
    num_threads : int
        Number of CPU threads.
    outfmt_fields : Sequence[str]
        BLAST output fields.

    Returns
    -------
    list[str]
        Command argument list suitable for subprocess execution.
    """
    cmd: List[str] = [
        "tblastn",
        "-task",
        task,
        "-query",
        str(query_fasta),
        "-db",
        str(db_path),
        "-out",
        str(out_tsv),
        "-evalue",
        str(evalue),
        "-word_size",
        str(word_size),
        "-seg",
        _normalize_seg(seg),
        "-soft_masking",
        _normalize_bool_flag(soft_masking),
        "-comp_based_stats",
        str(comp_based_stats),
        "-qcov_hsp_perc",
        str(qcov_hsp_perc),
        "-max_target_seqs",
        str(max_target_seqs),
        "-max_hsps",
        str(max_hsps),
        "-num_threads",
        str(num_threads),
        "-outfmt",
        _format_outfmt(outfmt_fields),
    ]

    if use_sw_tback:
        cmd.append("-use_sw_tback")

    return cmd


def _normalize_seg(seg: str | bool) -> str:
    """
    Normalize SEG filtering option for BLAST CLI usage.

    YAML configuration may automatically convert "yes"/"no"
    values into booleans. BLAST expects string literals
    ("yes" or "no") for the -seg option.

    Parameters
    ----------
    seg : str or bool
        SEG filtering configuration provided by pipeline input.

    Returns
    -------
    str
        BLAST-compatible SEG value.
    """

    if isinstance(seg, bool):
        return "yes" if seg else "no"

    return str(seg).lower()


def _normalize_bool_flag(value: bool) -> str:
    """
    Normalize boolean flag into BLAST CLI string representation.

    BLAST command-line options require lowercase string literals
    ("true" or "false") instead of Python boolean values.

    Parameters
    ----------
    value : bool
        Boolean flag value.

    Returns
    -------
    str
        BLAST-compatible boolean string.
    """

    return "true" if value else "false"
