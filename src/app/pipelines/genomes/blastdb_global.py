# pylint: disable=line-too-long, too-many-locals, broad-exception-caught, too-many-statements
"""
Global BLAST database pipeline.

This module orchestrates the end-to-end workflow for constructing
and registering a global BLAST database from a resolved genome scope.

Pipeline responsibilities:
1. Validate and summarize input scope
2. Build global BLAST database (processing stage)
3. Register BLAST database artifacts and provenance
4. Render execution summary
"""

from __future__ import annotations

# Standard Library Imports
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

# App Imports
from app.exceptions import PipelineError
from app.inputs.genus import load_genus_list
from app.processing.genomes.blastdb.build import build_global_blast_db
from app.processing.genomes.blastdb.global_scope import resolve_global_genome_scope
from app.registry.blastdb import is_blastdb_registered, load_blastdb, register_blastdb
from app.runtime.context import ExecutionContext
from app.runtime.logging.renderers import (
    render_batch_footer,
    render_batch_header,
    render_pipeline_error,
    render_stage_summary,
)
from app.runtime.logging.utils import format_duration, format_timestamp


# Public Pipeline API
def build_global_blastdb_pipeline(
    *,
    ctx: ExecutionContext,
    inputs: Dict,
    parameters: Dict,
    outputs: Dict,
    execution: Dict,
) -> Dict:
    """
    Build and register a global BLAST database.

    Parameters
    ----------
    ctx : ExecutionContext
        Runtime execution context for the current pipeline run.
    inputs : dict
        Pipeline input contract. Expected keys:
        - dataset_root : root directory containing registered
                         genus-level genome datasets
        - genera_list_path : path to genus list file
    parameters : dict
        Algorithm parameters controlling BLAST DB construction.
        Expected keys (optional):
        - db_name : str
        - molecule_type : {"nucl", "prot"}
    outputs : dict
        Pipeline output contract. Expected keys:
        - blastdb_root : root directory where BLAST DB artifacts
                         and registry files will be written
    execution : dict
        Execution policy controlling runtime behavior.
        Expected keys:
        - overwrite : bool

    Returns
    -------
    dict
        Pipeline execution summary.
    """
    # Initialize logger
    logger = ctx.logger

    # Parse execution parameters
    overwrite = execution.get("overwrite", False)

    # Parse inputs
    dataset_root = Path(inputs["dataset_root"])
    genera_list_path = Path(inputs["genera_list_path"])

    # Parse parameters
    db_name = parameters.get("db_name", "ampire_global")
    molecule_type = parameters.get("molecule_type", "nucl")

    # Parse outputs
    blastdb_root = Path(outputs["blastdb_root"])

    # Load genus list
    genus_list = load_genus_list(genera_list_path)

    # Prepare batch logging context
    pipeline_context = {
        "Input Summary": {
            "Dataset root": str(dataset_root),
            "Genus list path": str(genera_list_path),
            "Total genera requested": len(genus_list),
        },
        "BLAST DB Parameters": {
            "DB name": db_name,
            "Molecule type": molecule_type,
        },
        "Execution Policy": {
            "Overwrite enabled": overwrite,
        },
        "Output Summary": {
            "BLAST DB root": str(blastdb_root),
        },
    }

    started_at = datetime.now(timezone.utc)
    render_batch_header(
        logger=logger,
        title="Global BLAST DB Pipeline",
        context=pipeline_context,
    )

    logger.info("[ 'Global BLAST Database Construction' ]")

    # Track timing
    scope_started_at = None
    scope_finished_at = None
    build_started_at = None
    build_finished_at = None
    register_started_at = None
    register_finished_at = None

    try:
        # Stage 0: Check registry reuse shortcut
        reuse_started_at = datetime.now(timezone.utc)

        if (
            is_blastdb_registered(db_name=db_name, data_root=blastdb_root)
            and not overwrite
        ):
            blastdb_summary = load_blastdb(db_name=db_name, data_root=blastdb_root)

            logger.add_divider(length=120, border="+", fill="-")
            render_stage_summary(
                logger=logger,
                title="[ Reuse BLAST DB ]",
                summary_lines=[
                    f"▸ DB name       : '{db_name}'",
                    f"▸ Molecule type : '{blastdb_summary.get('blastdb', {}).get('type', 'unknown')}'",
                    "▸ Build state   : 'reused'",
                    "▸ Registry state: 'reused'",
                ],
                detail_lines=[
                    "  - Existing registry artifact detected",
                    "  - Overwrite is disabled",
                    f"      • Registry YAML : '{blastdb_root / db_name / 'blastdb.yaml'}'",
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
                title="Global BLAST DB Pipeline",
                context=summary_context,
            )

            return {
                "status": status,
                "blastdb": db_name,
                "blastdb_dir": str(blastdb_root / db_name),
                "metadata": str(blastdb_root / db_name / "blastdb.yaml"),
                "duration": format_duration(total_duration),
            }

        # Stage 1: Resolve global genome scope
        scope_started_at = datetime.now(timezone.utc)

        scope = resolve_global_genome_scope(
            dataset_root=dataset_root,
            genus_list=genus_list,
        )

        excluded = scope.get("excluded", [])
        summary = scope.get("summary", {})

        logger.add_divider(length=120, border="+", fill="-")
        render_stage_summary(
            logger=logger,
            title="[ Resolve global genome scope ]",
            summary_lines=[
                f"▸ Total genera requested : {summary.get('total_genera')}",
                f"▸ Genera included        : {summary.get('included_genera')}",
                f"▸ Genera excluded        : {summary.get('excluded_genera')}",
                f"▸ Total genomes included : {summary.get('total_genomes')}",
            ],
            detail_lines=[
                "  - Excluded genera :" if excluded else "  - No excluded genera",
                *[f"      • {rec['genus']} ('{rec['reason']}')" for rec in excluded],
            ],
        )

        scope_finished_at = datetime.now(timezone.utc)

        # Stage 2: Build global BLAST DB
        build_started_at = datetime.now(timezone.utc)

        blastdb_summary = build_global_blast_db(
            scope=scope,
            output_dir=blastdb_root,
            db_name=db_name,
            molecule_type=molecule_type,
            overwrite=overwrite,
        )

        build_state = blastdb_summary.get("blastdb", {}).get("status", "unknown")

        logger.add_divider(length=120, border="+", fill="-")
        render_stage_summary(
            logger=logger,
            title="[ Build BLAST DB ]",
            summary_lines=[
                f"▸ DB name       : '{db_name}'",
                f"▸ Molecule type : '{molecule_type}'",
                f"▸ Build state   : '{build_state}'",
            ],
            detail_lines=[
                "  - Artifacts :",
                f"      • FASTA : '{blastdb_summary['artifacts']['fasta_path']}'",
                f"      • DB    : '{blastdb_summary['artifacts']['db_base']}'",
            ],
        )

        build_finished_at = datetime.now(timezone.utc)

        # Stage 3: Register BLAST DB
        register_started_at = datetime.now(timezone.utc)

        blastdb_yaml = register_blastdb(
            db_name=db_name,
            data_root=blastdb_root,
            blastdb_summary=blastdb_summary,
            overwrite=overwrite,
        )

        registry_state = "registered"
        registry_details = [
            "  - Artifacts :",
            f"      • Registry YAML : '{blastdb_yaml}'",
        ]

        logger.add_divider(length=120, border="+", fill="-")
        render_stage_summary(
            logger=logger,
            title="[ Register BLAST DB ]",
            summary_lines=[
                f"▸ Registry state : '{registry_state}'",
            ],
            detail_lines=registry_details,
        )

        register_finished_at = datetime.now(timezone.utc)

        # Final timing summary
        logger.add_divider(length=120, border="+", fill="-")
        timing_lines = [
            "[ Timing Summary ]",
            f"▸ Resolve scope   : {format_duration(scope_finished_at - scope_started_at)}",
            f"▸ Build BLAST DB  : {format_duration(build_finished_at - build_started_at)}",
            f"▸ Register DB     : {format_duration(register_finished_at - register_started_at)}",
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

    except PipelineError as exc:
        render_pipeline_error(logger=logger, error=exc)
        status = "failed"
        reason = exc.reason

    except Exception as exc:
        fallback_error = PipelineError(
            pipeline="genomes.blastdb_global",
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
        title="Global BLAST DB Pipeline",
        context=summary_context,
    )

    return {
        "status": status,
        "blastdb": db_name,
        "blastdb_dir": str(blastdb_root / db_name),
        "metadata": str(blastdb_root / db_name / "blastdb.yaml"),
        "duration": format_duration(total_duration),
    }
