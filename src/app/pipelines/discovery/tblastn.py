# pylint: disable=line-too-long, too-many-locals, broad-exception-caught, too-many-statements
"""
TBLASTN discovery search pipeline.

This pipeline orchestrates genome-scale AMP discovery searches using
tblastn against a pre-built nucleotide BLAST database.

Stages:
1. Execute tblastn search
2. Normalize alignment output
3. Compute search fingerprint
4. Register discovery dataset (registry stage)
"""

from __future__ import annotations

# Standard Library
import logging
import shutil
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

# App Imports
from app.exceptions import PipelineError
from app.processing.discovery.tblastn.fingerprint import (
    SCIENTIFIC_PARAMETER_KEYS,
    compute_search_fingerprint,
    compute_search_id,
)
from app.processing.discovery.tblastn.normalize import (
    CANONICAL_COLUMNS,
    normalize_tblastn_output,
)
from app.processing.discovery.tblastn.run_search import run_tblastn_search
from app.registry.discovery import (
    is_discovery_registered,
    load_discovery,
    register_discovery,
)
from app.runtime.context import ExecutionContext
from app.runtime.logging.renderers import (
    render_batch_footer,
    render_batch_header,
    render_pipeline_error,
    render_stage_summary,
)
from app.runtime.logging.utils import format_duration, format_timestamp


# Public Pipeline API
def tblastn_discovery_search(
    *,
    ctx: ExecutionContext,
    inputs: Dict,
    parameters: Dict,
    outputs: Dict,
    execution: Dict,
) -> Dict:
    """
    Execute a tblastn discovery search pipeline.

    Returns
    -------
    Dict
        Execution summary.
    """

    # Initialize logger
    logger = ctx.logger

    # Parse execution parameters
    overwrite = execution.get("overwrite", False)

    # Parse inputs
    query_fasta = Path(inputs["query_fasta"])
    db_name = inputs["db_name"]
    blastdb_root = Path(inputs["blastdb_root"])
    db_prefix = blastdb_root / db_name / db_name

    # Parse parameters
    search_series = parameters.get("search_series")
    search_label = parameters.get("search_label")
    evalue = parameters.get("evalue", 10.0)
    word_size = parameters.get("word_size", 3)
    seg = parameters.get("seg", "yes")
    soft_masking = parameters.get("soft_masking", True)
    comp_based_stats = parameters.get("comp_based_stats", 2)
    use_sw_tback = parameters.get("use_sw_tback", False)
    max_target_seqs = parameters.get("max_target_seqs", 500)
    max_hsps = parameters.get("max_hsps", 1)
    qcov_hsp_perc = parameters.get("qcov_hsp_perc", 0)
    batch_size = parameters.get("batch_size", 100)
    num_threads = parameters.get("num_threads", 10)  # runtime-only

    # Parse outputs
    discovery_root = Path(outputs["discovery_root"])
    discovery_root.mkdir(parents=True, exist_ok=True)

    # Prepare logging context
    pipeline_context = {
        "Input Summary": {
            "Query FASTA": str(query_fasta),
            "BLAST DB name": db_name,
            "BLAST DB prefix": str(db_prefix),
        },
        "Search Identity Parameters": {
            "Search series": search_series,
            "Search label": search_label,
            "Batch size": batch_size,
            "E-value": evalue,
            "Word size": word_size,
            "SEG filtering": seg,
            "Soft masking": soft_masking,
            "Composition stats": comp_based_stats,
            "Smith-Waterman": use_sw_tback,
            "Max target seqs": max_target_seqs,
            "Max HSPs": max_hsps,
            "Query coverage": qcov_hsp_perc,
            "Threads": num_threads,
        },
        "Execution Policy": {
            "Overwrite enabled": overwrite,
        },
        "Output Summary": {
            "Discovery root": str(discovery_root),
        },
    }

    started_at = datetime.now(timezone.utc)
    render_batch_header(
        logger=logger,
        title="TBLASTN Discovery Search Pipeline",
        context=pipeline_context,
    )

    logger.info("[ 'TBLASTN Discovery Search' ]")

    search_id = None
    discovery_yaml = None

    # Track timing
    search_started_at = None
    search_finished_at = None
    normalize_started_at = None
    normalize_finished_at = None
    register_started_at = None
    register_finished_at = None

    try:
        # Stage 0: Validate + canonicalize identity parameters
        identity_params = _validate_and_canonicalize_identity_parameters(parameters)

        # Defensive merge: keep original config but ensure identity keys are canonical
        parameters = {**parameters, **identity_params}

        # Stage 0: Check registry reuse shortcut
        reuse_started_at = datetime.now(timezone.utc)

        fingerprint_preview = compute_search_fingerprint(
            query_fasta=query_fasta,
            db_prefix=db_prefix,
            outfmt_fields=CANONICAL_COLUMNS,
            parameters=identity_params,
        )

        search_id_preview = compute_search_id(fingerprint_preview)

        # Create isolated execution workspace
        run_uuid = uuid.uuid4().hex[:8]
        tmp_dir = discovery_root / "_tmp" / search_id_preview / f"run_{run_uuid}"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        raw_output = tmp_dir / "alignments_raw.tsv"
        normalized_output = tmp_dir / "alignments.tsv"

        if (
            is_discovery_registered(
                search_id=search_id_preview,
                data_root=discovery_root,
            )
            and not overwrite
        ):

            discovery_summary = load_discovery(
                search_id=search_id_preview,
                data_root=discovery_root,
            )

            logger.add_divider(length=120, border="+", fill="-")
            render_stage_summary(
                logger=logger,
                title="[ Reuse Discovery Dataset ]",
                summary_lines=[
                    f"▸ Search ID     : '{search_id_preview}'",
                    "▸ Execution     : 'reused'",
                    "▸ Registry      : 'reused'",
                ],
                detail_lines=[
                    "  - Existing discovery registry detected",
                    "  - Overwrite is disabled",
                    f"      • Registry YAML : '{discovery_root / search_id_preview / 'discovery.yaml'}'",
                ],
            )

            reuse_finished_at = datetime.now(timezone.utc)

            # Timing summary (reuse mode)
            logger.add_divider(length=120, border="+", fill="-")
            timing_lines = [
                "[ Timing Summary ]",
                f"▸ Reuse registry : {format_duration(reuse_finished_at - reuse_started_at)}",
            ]

            logger.log_with_borders(
                level=logging.INFO,
                message="\n".join(timing_lines),
                length=120,
            )
            logger.add_divider(length=120, border="+", fill="-")
            logger.add_spacer(level=logging.INFO, lines=1)

            status = "success"
            reason = None

            finished_at = datetime.now(timezone.utc)
            total_duration = finished_at - started_at

            summary_context = {
                "Results": {
                    "Status": status,
                    "Reason": reason,
                },
                "Timing": {
                    "Start Time": format_timestamp(started_at),
                    "End Time": format_timestamp(finished_at),
                    "Duration": format_duration(total_duration),
                },
            }

            render_batch_footer(
                logger=logger,
                title="TBLASTN Discovery Search Pipeline",
                context=summary_context,
            )

            return {
                "status": status,
                "search_id": search_id_preview,
                "metadata": str(discovery_root / search_id_preview / "discovery.yaml"),
            }

        # Stage 1: Run tblastn (batched)
        search_started_at = datetime.now(timezone.utc)

        batch_dir = tmp_dir / "query_batches"
        batch_dir.mkdir(parents=True, exist_ok=True)

        batch_files = _split_fasta(query_fasta, batch_dir, batch_size)
        total_batches = len(batch_files)
        if total_batches == 0:
            raise PipelineError(
                pipeline="discovery.tblastn",
                stage="split_query_fasta",
                reason="Query FASTA produced zero batches.",
                action="Check FASTA formatting or batch_size parameter.",
            )

        logger.add_divider(length=120, border="+", fill="-")
        render_stage_summary(
            logger=logger,
            title="[ Split Query FASTA ]",
            summary_lines=[
                f"▸ Query FASTA   : '{query_fasta.name}'",
                f"▸ Batch size    : {batch_size}",
                f"▸ Total batches : {total_batches}",
            ],
            detail_lines=[
                f"  - Batch directory : '{batch_dir}'",
            ],
        )
        logger.add_divider(length=120, border="+", fill="-")

        raw_parts = []

        for i, batch_file in enumerate(batch_files, start=1):

            batch_output = tmp_dir / f"alignments_raw_{i-1:04d}.tsv"

            batch_start = datetime.now(timezone.utc)

            run_tblastn_search(
                query_fasta=batch_file,
                db_path=db_prefix,
                out_tsv=batch_output,
                evalue=identity_params["evalue"],
                word_size=identity_params["word_size"],
                seg=identity_params["seg"],
                soft_masking=identity_params["soft_masking"],
                comp_based_stats=identity_params["comp_based_stats"],
                use_sw_tback=identity_params["use_sw_tback"],
                max_target_seqs=identity_params["max_target_seqs"],
                max_hsps=identity_params["max_hsps"],
                qcov_hsp_perc=identity_params["qcov_hsp_perc"],
                num_threads=num_threads,
                overwrite=True,
            )

            batch_end = datetime.now(timezone.utc)
            logger.log_with_borders(
                level=logging.INFO,
                message="\n".join(
                    [
                        f"[ Batch {i}/{total_batches} ] Running tblastn on '{batch_file.name}'. "
                        f"Completed in {format_duration(batch_end - batch_start)}",
                    ]
                ),
                length=120,
            )

            raw_parts.append(batch_output)

        # Merge outputs
        with open(raw_output, "w", encoding="utf-8") as outfile:
            for part in raw_parts:
                with open(part, "r", encoding="utf-8") as infile:
                    shutil.copyfileobj(infile, outfile)

        search_finished_at = datetime.now(timezone.utc)

        # Stage 2: Normalize output
        normalize_started_at = datetime.now(timezone.utc)

        df = normalize_tblastn_output(raw_output)
        df.to_csv(normalized_output, sep="\t", index=False)

        normalize_finished_at = datetime.now(timezone.utc)

        # Stage 3: Fingerprint
        fingerprint = compute_search_fingerprint(
            query_fasta=query_fasta,
            db_prefix=db_prefix,
            outfmt_fields=CANONICAL_COLUMNS,
            parameters=identity_params,
        )

        search_id = compute_search_id(fingerprint)

        # Final dataset directory
        final_dir = discovery_root / search_id
        final_dir.mkdir(parents=True, exist_ok=True)

        total_hsps = len(df)
        total_queries = df["qseqid"].nunique() if not df.empty else 0

        # Move artifacts into dataset directory
        final_raw = final_dir / "alignments_raw.tsv"
        final_norm = final_dir / "alignments.tsv"

        # Commit guard: another execution may have finished first
        if final_raw.exists() and not overwrite:

            registry_state = "reused"

            logger.add_divider(length=120, border="+", fill="-")
            render_stage_summary(
                logger=logger,
                title="[ Reuse Discovery Results (Concurrent Commit) ]",
                summary_lines=[
                    f"▸ Search ID      : '{search_id}'",
                    "▸ Registry state : 'reused'",
                ],
                detail_lines=[
                    "  - Another execution completed first",
                    "  - Current run reused finalized dataset",
                ],
            )
            logger.add_divider(length=120, border="+", fill="-")
            logger.add_spacer(level=logging.INFO, lines=1)

            discovery_yaml = final_dir / "discovery.yaml"

            status = "success"
            reason = None

            _cleanup_tmp_workspace(tmp_dir)

            # Skip remaining stages but KEEP footer
            raise StopIteration

        raw_output.replace(final_raw)
        normalized_output.replace(final_norm)

        raw_output = final_raw
        normalized_output = final_norm

        # Stage 4: Register discovery results
        register_started_at = datetime.now(timezone.utc)
        discovery_summary = {
            "discovery": {
                "search_id": search_id,
                "status": "completed",
                "created_at": format_timestamp(datetime.now(timezone.utc)),
            },
            "fingerprint": fingerprint,
            "inputs": {
                "query_fasta": str(query_fasta),
                "blastdb_prefix": str(db_prefix),
            },
            "parameters_raw": {
                k: parameters.get(k)
                for k in SCIENTIFIC_PARAMETER_KEYS
                if k in parameters
            },
            "parameters": {
                k: identity_params.get(k)
                for k in SCIENTIFIC_PARAMETER_KEYS
                if k in identity_params
            },
            "statistics": {
                "total_queries": total_queries,
                "total_hsps": total_hsps,
            },
            "artifacts": {
                "raw_alignments": str(raw_output),
                "normalized_alignments": str(normalized_output),
            },
        }
        discovery_yaml = register_discovery(
            search_id=search_id,
            data_root=discovery_root,
            discovery_summary=discovery_summary,
            overwrite=overwrite,
        )

        registry_state = "registered"
        registry_details = [
            "  - Artifacts :",
            f"      • Registry YAML : '{discovery_yaml}'",
        ]

        logger.add_divider(length=120, border="+", fill="-")
        render_stage_summary(
            logger=logger,
            title="[ Register Discovery Results ]",
            summary_lines=[
                f"▸ Search ID      : '{search_id}'",
                f"▸ Registry state : '{registry_state}'",
            ],
            detail_lines=registry_details,
        )

        register_finished_at = datetime.now(timezone.utc)

        # Cleanup temporary workspace after successful commit
        _cleanup_tmp_workspace(tmp_dir)

        # Final timing summary
        logger.add_divider(length=120, border="+", fill="-")
        timing_lines = [
            "[ Timing Summary ]",
            f"▸ tblastn search : {format_duration(search_finished_at - search_started_at)}",
            f"▸ normalization  : {format_duration(normalize_finished_at - normalize_started_at)}",
            f"▸ register       : {format_duration(register_finished_at - register_started_at)}",
        ]

        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(timing_lines),
            length=120,
        )
        logger.add_divider(length=120, border="+", fill="-")
        logger.add_spacer(level=logging.INFO, lines=1)

        status = "success"
        reason = None

    except StopIteration:
        pass

    except PipelineError as exc:
        render_pipeline_error(logger=logger, error=exc)
        status = "failed"
        reason = exc.reason

    except Exception as exc:
        fallback_error = PipelineError(
            pipeline="discovery.tblastn",
            stage="unknown",
            reason=str(exc),
            action="Inspect stack trace and pipeline configuration.",
        )
        render_pipeline_error(logger=logger, error=fallback_error)
        status = "failed"
        reason = str(exc)

    # Pipeline summary
    finished_at = datetime.now(timezone.utc)
    total_duration = finished_at - started_at

    summary_context = {
        "Results": {"Status": status, "Reason": reason},
        "Timing": {
            "Start Time": format_timestamp(started_at),
            "End Time": format_timestamp(finished_at),
            "Duration": format_duration(total_duration),
        },
    }

    render_batch_footer(
        logger=logger,
        title="TBLASTN Discovery Search Pipeline",
        context=summary_context,
    )

    return {
        "status": status,
        "search_id": search_id if status == "success" else None,
        "metadata": str(discovery_yaml) if status == "success" else None,
        "duration": format_duration(total_duration),
    }


# Internal Helpers
def _cleanup_tmp_workspace(tmp_dir: Path) -> None:
    """
    Remove temporary execution workspace if empty.

    This function safely removes:
    - run_<uuid> directory
    - parent <search_id> directory (if empty)

    The global "_tmp" directory is never removed.

    Parameters
    ----------
    tmp_dir : pathlib.Path
        Temporary run workspace directory.
    """

    try:
        # Remove run directory
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)

        # Remove search_id directory if empty
        search_tmp_dir = tmp_dir.parent

        if search_tmp_dir.exists() and not any(search_tmp_dir.iterdir()):
            search_tmp_dir.rmdir()

    except OSError:
        # Directory not empty or concurrent run still active
        pass


def _validate_and_canonicalize_identity_parameters(parameters: Dict) -> Dict:
    """
    Validate and canonicalize tblastn search identity parameters.

    This function enforces deterministic parameter types and representations
    for fingerprint/search_id stability.

    Notes
    -----
    This function only targets identity parameters used for fingerprinting.
    Runtime-only fields (e.g., num_threads) are intentionally excluded.

    Parameters
    ----------
    parameters : Dict
        Raw parameter mapping parsed from YAML/CLI.

    Returns
    -------
    Dict
        Canonicalized identity parameter mapping.

    Raises
    ------
    PipelineError
        If any identity parameter is invalid or cannot be canonicalized.
    """

    try:
        canonical: Dict[str, Any] = {}

        # Scientific identity parameters
        canonical["evalue"] = _coerce_float(
            parameters.get("evalue", 10.0), name="evalue"
        )
        canonical["word_size"] = _coerce_int(
            parameters.get("word_size", 3), name="word_size"
        )
        canonical["seg"] = _coerce_seg(parameters.get("seg", "yes"))
        canonical["soft_masking"] = _coerce_bool(
            parameters.get("soft_masking", True), name="soft_masking"
        )
        canonical["comp_based_stats"] = _coerce_int(
            parameters.get("comp_based_stats", 2), name="comp_based_stats"
        )
        canonical["use_sw_tback"] = _coerce_bool(
            parameters.get("use_sw_tback", False), name="use_sw_tback"
        )
        canonical["max_target_seqs"] = _coerce_int(
            parameters.get("max_target_seqs", 500), name="max_target_seqs"
        )
        canonical["max_hsps"] = _coerce_int(
            parameters.get("max_hsps", 1), name="max_hsps"
        )
        canonical["qcov_hsp_perc"] = _coerce_float(
            parameters.get("qcov_hsp_perc", 0),
            name="qcov_hsp_perc",
        )

        # Optional identity labels (kept stable as strings)
        if parameters.get("search_series") is not None:
            canonical["search_series"] = str(parameters.get("search_series") or "")
        if parameters.get("search_label") is not None:
            canonical["search_label"] = str(parameters.get("search_label") or "")

        # Identity-level validation rules
        if canonical["evalue"] <= 0:
            raise ValueError("evalue must be > 0.")
        if canonical["word_size"] <= 0:
            raise ValueError("word_size must be > 0.")
        if canonical["comp_based_stats"] < 0:
            raise ValueError("comp_based_stats must be >= 0.")
        if canonical["max_target_seqs"] <= 0:
            raise ValueError("max_target_seqs must be > 0.")
        if canonical["max_hsps"] <= 0:
            raise ValueError("max_hsps must be > 0.")
        if canonical["qcov_hsp_perc"] < 0:
            raise ValueError("qcov_hsp_perc must be >= 0.")

        return canonical

    except Exception as exc:
        raise PipelineError(
            pipeline="discovery.tblastn",
            stage="validate_identity_parameters",
            reason=str(exc),
            action="Fix parameter types/values in the YAML configuration.",
        ) from exc


def _coerce_float(value: Any, *, name: str) -> float:
    """
    Coerce numeric parameter into float.

    Parameters
    ----------
    value : Any
        Raw value parsed from YAML/CLI.
    name : str
        Parameter name used for error messages.

    Returns
    -------
    float
        Coerced float value.

    Raises
    ------
    ValueError
        If value cannot be coerced into float.
    """
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a float, got bool.")
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        return float(value.strip())
    raise ValueError(f"{name} must be a float, got {type(value).__name__}.")


def _coerce_int(value: Any, *, name: str) -> int:
    """
    Coerce numeric parameter into int.

    Parameters
    ----------
    value : Any
        Raw value parsed from YAML/CLI.
    name : str
        Parameter name used for error messages.

    Returns
    -------
    int
        Coerced integer value.

    Raises
    ------
    ValueError
        If value cannot be coerced into int.
    """
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an int, got bool.")
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not value.is_integer():
            raise ValueError(f"{name} must be an int, got float '{value}'.")
        return int(value)
    if isinstance(value, str):
        s = value.strip()
        try:
            return int(s)
        except ValueError:
            f = float(s)
            if not f.is_integer():
                raise ValueError(f"{name} must be an int, got '{value}'.") from None
            return int(f)
    raise ValueError(f"{name} must be an int, got {type(value).__name__}.")


def _coerce_bool(value: Any, *, name: str) -> bool:
    """
    Coerce YAML/CLI value into bool.

    Parameters
    ----------
    value : Any
        Raw value parsed from YAML/CLI.
    name : str
        Parameter name used for error messages.

    Returns
    -------
    bool
        Coerced boolean value.

    Raises
    ------
    ValueError
        If value cannot be coerced into bool.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        if value in (0, 1):
            return bool(value)
        raise ValueError(f"{name} must be a bool, got int '{value}'.")
    if isinstance(value, str):
        s = value.strip().lower()
        if s in ("true", "yes", "y", "1"):
            return True
        if s in ("false", "no", "n", "0"):
            return False
        raise ValueError(f"{name} must be a bool, got '{value}'.")
    raise ValueError(f"{name} must be a bool, got {type(value).__name__}.")


def _coerce_seg(seg: Any) -> str:
    """
    Coerce SEG filtering option into BLAST-compatible string.

    Parameters
    ----------
    seg : Any
        Raw seg value parsed from YAML/CLI.

    Returns
    -------
    str
        Normalized seg string ("yes" or "no").

    Raises
    ------
    ValueError
        If seg cannot be coerced into a valid value.
    """
    if isinstance(seg, bool):
        return "yes" if seg else "no"

    if isinstance(seg, str):
        s = seg.strip().lower()
        if s in ("yes", "true", "y", "1"):
            return "yes"
        if s in ("no", "false", "n", "0"):
            return "no"
        raise ValueError(f"seg must be 'yes' or 'no', got '{seg}'.")

    raise ValueError(f"seg must be 'yes' or 'no', got {type(seg).__name__}.")


def _split_fasta(
    fasta_path: Path,
    output_dir: Path,
    batch_size: int,
) -> List[Path]:
    """
    Split a FASTA file into multiple smaller batch files.

    Parameters
    ----------
    fasta_path : Path
        Input FASTA file containing query sequences.
    output_dir : Path
        Directory where batch FASTA files will be written.
    batch_size : int
        Maximum number of sequences per batch.

    Returns
    -------
    list[Path]
        List of generated batch FASTA file paths.
    """

    batches: List[Path] = []
    batch_index = 0
    seq_count = 0
    current_batch: List[Tuple[str, str]] = []

    header: str | None = None
    seq_lines: List[str] = []

    with open(fasta_path, "r", encoding="utf-8") as f:

        for line in f:
            line = line.rstrip()

            if line.startswith(">"):

                # Flush previous record
                if header is not None:
                    current_batch.append((header, "".join(seq_lines)))
                    seq_count += 1

                    # Emit batch if full
                    if seq_count == batch_size:
                        batch_file = output_dir / f"batch_{batch_index:04d}.fasta"
                        _write_batch(batch_file, current_batch)
                        batches.append(batch_file)

                        batch_index += 1
                        current_batch = []
                        seq_count = 0

                header = line
                seq_lines = []

            else:
                seq_lines.append(line)

        # Flush last record
        if header is not None:
            current_batch.append((header, "".join(seq_lines)))
            seq_count += 1

    # Write remaining records
    if current_batch:
        batch_file = output_dir / f"batch_{batch_index:04d}.fasta"
        _write_batch(batch_file, current_batch)
        batches.append(batch_file)

    return batches


def _write_batch(
    path: Path,
    records: List[Tuple[str, str]],
) -> None:
    """
    Write a list of FASTA records to file.

    Parameters
    ----------
    path : Path
        Output FASTA file path.
    records : list[tuple[str, str]]
        FASTA records as (header, sequence) tuples.

    Returns
    -------
    None
    """

    with open(path, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f"{header}\n{seq}\n")
