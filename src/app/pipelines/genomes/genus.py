# pylint: disable=too-many-locals, broad-exception-caught, too-many-arguments, too-many-statements
"""
Genus-level genome dataset pipeline.

This module orchestrates the end-to-end workflow for building a
reproducible genome dataset for a given genus, including:

1. Ensuring taxonomy is registered
2. Resolving species TaxIDs
3. Downloading genome assemblies
4. Registering dataset artifacts and provenance
"""

from __future__ import annotations

# Standard Library Imports
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

# App Imports
from app.exceptions import PipelineError
from app.inputs.genus import load_genus_list
from app.processing.genomes.download.by_genus import download_genomes_from_genus
from app.processing.genomes.taxonomy.species_from_genus import (
    resolve_species_from_genus,
)
from app.registry.genomes import (
    is_genus_genomes_registered,
    register_genus_genomes_dataset,
)
from app.registry.taxonomy import register_genus_taxonomy
from app.runtime.context import ExecutionContext
from app.runtime.logging.renderers import (
    render_batch_footer,
    render_batch_header,
    render_pipeline_error,
    render_stage_summary,
)
from app.runtime.logging.utils import format_duration, format_timestamp


# Public Pipeline API
def build_genus_genomes_batch_dataset(
    *,
    ctx: ExecutionContext,
    inputs: Dict,
    parameters: Dict,
    outputs: Dict,
    execution: Dict,
) -> Dict:
    """
    Build genus-level genome datasets in batch mode.

    This pipeline orchestrates the execution of the single-genus
    genome dataset pipeline across multiple genera.

    The batch pipeline is responsible for:
    - Loading genus list
    - Iterative pipeline execution
    - Error isolation and handling
    - Execution summary reporting

    Parameters
    ----------
    ctx : ExecutionContext
        Runtime execution context for the current pipeline run.
    inputs : dict
        Pipeline input contract. Expected keys:
        - genus_list_path : path to genus list file
    parameters : dict
        Algorithm and retrieval parameters controlling scientific
        behavior of the pipeline.
    outputs : dict
        Pipeline output contract. Expected keys:
        - root : root directory where all dataset artifacts
                 produced by this pipeline will be written.
    execution : dict
        Execution policy controlling runtime behavior.
        Expected keys:
        - overwrite : bool

    Returns
    -------
    Dict
        Batch execution summary with the following structure:
        {
            "total": int,
            "success": int,
            "failed": int,
            "skipped": int,
            "results": [
                {
                    "genus": str,
                    "status": "success" | "failed" | "skipped",
                    "stage": Optional[str],
                    "reason": Optional[str],
                    "genomes_root": Optional[str],
                },
                ...
            ]
        }
    """
    # Initialize logger
    logger = ctx.logger

    # Parse execution parameters
    overwrite = execution.get("overwrite", False)
    genus_list_path = Path(inputs["genus_list_path"])
    dataset_root = Path(outputs["root"])

    # Parse algorithm parameters
    assembly_level = parameters.get("assembly_level", "complete,chromosome")
    refseq_categories = parameters.get("refseq_categories", "all")
    source = parameters.get("source", "refseq")
    parallel = parameters.get("parallel", 10)
    batch_size = parameters.get("batch_size", 3000)

    # Load genus list
    genus_list = load_genus_list(genus_list_path)
    total = len(genus_list)

    # Prepare batch logging context
    pipeline_context = {
        "Input Summary": {
            "Genus list path": str(genus_list_path),
            "Total genera": total,
        },
        "Retrieval Parameters": {
            "Source": source,
            "Assembly level": assembly_level,
            "RefSeq categories": refseq_categories,
            "Parallel workers": parallel,
            "Batch size": batch_size,
        },
        "Execution Policy": {
            "Overwrite enabled": overwrite,
        },
        "Output Summary": {
            "Dataset root": str(dataset_root),
        },
    }

    # Record batch start time
    batch_started_at = datetime.now(timezone.utc)

    # Log batch initialization details
    render_batch_header(
        logger=logger,
        title="Genus Genome Batch Pipeline",
        context=pipeline_context,
    )

    # Step 1. Initialize batch result container
    results = []
    success = 0
    failed = 0
    skipped = 0

    # Step 2. Iterate through genera
    for index, genus_name in enumerate(genus_list, start=1):

        # Record genus start time
        genus_started_at = datetime.now(timezone.utc)

        # Log current genus processing
        logger.info(f"[ {index} / {total} ] Genus: '{genus_name}'")
        logger.add_divider(length=120, border="+", fill="-")
        try:
            # Build genus-level genome dataset
            genomes_root = _build_genus_genomes_dataset(
                ctx=ctx,
                genus_name=genus_name,
                dataset_root=dataset_root,
                assembly_level=assembly_level,
                refseq_categories=refseq_categories,
                source=source,
                parallel=parallel,
                batch_size=batch_size,
                overwrite=overwrite,
            )

            # Log success and record result
            results.append(
                {
                    "genus": genus_name,
                    "status": "success",
                    "genomes_root": str(genomes_root),
                    "stage": None,
                    "reason": None,
                }
            )
            success += 1

        except PipelineError as exc:

            # Log failure and record result
            render_pipeline_error(
                logger=logger,
                error=exc,
            )
            results.append(
                {
                    "genus": genus_name,
                    "status": "failed",
                    "stage": exc.stage,
                    "reason": exc.reason,
                    "genomes_root": None,
                }
            )
            failed += 1

        except Exception as exc:

            # Log unexpected error and record result
            fallback_error = PipelineError(
                pipeline="genomes.genus",
                stage="unknown",
                reason=str(exc),
                action="Inspect stack trace and pipeline configuration.",
            )
            render_pipeline_error(
                logger=logger,
                error=fallback_error,
            )
            results.append(
                {
                    "genus": genus_name,
                    "status": "failed",
                    "stage": "unknown",
                    "reason": str(exc),
                    "action": fallback_error.action,
                    "genomes_root": None,
                }
            )
            failed += 1

        finally:
            genus_finished_at = datetime.now(timezone.utc)
            genus_duration = genus_finished_at - genus_started_at
            logger.add_divider(length=120, border="+", fill="-")
            logger.log_with_borders(
                level=logging.INFO,
                message="[ Timing ] " f"Duration : {format_duration(genus_duration)}",
                length=120,
            )

        logger.add_divider(length=120, border="+", fill="-")
        logger.add_spacer(level=logging.INFO, lines=1)

    # Pipeline summary
    batch_finished_at = datetime.now(timezone.utc)
    duration = batch_finished_at - batch_started_at
    summary_context = {
        "Results": {
            "Total": total,
            "Success": success,
            "Failed": failed,
            "Skipped": skipped,
        },
        "Timing": {
            "Start Time": format_timestamp(batch_started_at),
            "End Time": format_timestamp(batch_finished_at),
            "Duration": format_duration(duration),
        },
    }
    render_batch_footer(
        logger=logger,
        title="Genus Genome Batch Pipeline",
        context=summary_context,
    )

    # Step 3. Return batch summary
    return {
        "total": len(genus_list),
        "success": success,
        "failed": failed,
        "skipped": skipped,
        "results": results,
    }


# Internal Helpers
def _build_genus_genomes_dataset(
    *,
    ctx: ExecutionContext,
    genus_name: str,
    dataset_root: Path,
    assembly_level: str,
    refseq_categories: str,
    source: str,
    parallel: int,
    batch_size: int,
    overwrite: bool = False,
) -> Optional[Path]:
    """
    Build a genus-level genome dataset.

    Parameters
    ----------
    ctx : ExecutionContext
        Runtime execution context inherited from the batch pipeline.
    genus_name : str
        Scientific genus name.
    dataset_root : pathlib.Path
        Root directory under which all genus-level dataset artifacts
        will be created.
    assembly_level : str
        Genome assembly levels to retrieve.
    refseq_categories : str
        RefSeq categories to include.
    source : str
        Genome source database (e.g., "refseq").
    parallel : int
        Number of parallel workers for genome download.
    batch_size : int
        Batch size for genome download.
    overwrite : bool, optional
        Whether to overwrite existing registered genome datasets.

    Returns
    -------
    Optional[pathlib.Path]
        Path to the registered genus-level genomes dataset directory.

    Raises
    ------
    PipelineError
        If any pipeline stage fails.
    """
    # Initialize logger
    logger = ctx.logger

    # Step 1: Resolve taxonomy (domain stage)
    genus_metadata, species_table = resolve_species_from_genus(genus_name)
    genus_taxid = genus_metadata["genus"].get("taxid")

    # Validate taxonomy resolution
    if genus_taxid is None:
        logger.log_with_borders(
            level=logging.INFO, message="[ Resolve taxonomy ]", length=120
        )
        raise PipelineError(
            pipeline="genomes.genus",
            stage="taxonomy",
            reason="NCBI taxonomy not found.",
            action="Please check genus spelling or taxonomy source.",
            context={"genus": genus_name},
        )

    if species_table.empty:
        logger.log_with_borders(
            level=logging.INFO, message="[ Resolve taxonomy ]", length=120
        )
        raise PipelineError(
            pipeline="genomes.genus",
            stage="taxonomy",
            reason="No species-level taxa found under this genus.",
            action="This genus exists but has no species-level descendants.",
            context={"genus": genus_name},
        )

    # Step 2: Register taxonomy (registry stage)
    try:
        genus_yaml, species_csv = register_genus_taxonomy(
            genus_name=genus_name,
            data_root=dataset_root,
            genus_metadata=genus_metadata,
            species_table=species_table,
            overwrite=overwrite,
        )
        taxonomy_state = "taxonomy_registered"
        taxonomy_details = [
            "  - Artifacts :",
            f"     • Genus metadata : '{genus_yaml}'",
            f"     • Species table  : '{species_csv}'",
        ]

    except FileExistsError as exc:
        taxonomy_state = "taxonomy_already_registered"
        taxonomy_details = [
            f"  - {str(exc)}",
        ]

    # Log taxonomy resolution summary
    render_stage_summary(
        logger=logger,
        title="[ Resolve taxonomy ]",
        summary_lines=[
            f"▸ Genus taxid        : {genus_taxid}",
            f"▸ Species discovered : {len(species_table)}",
            f"▸ Registry state     : '{taxonomy_state}'",
        ],
        detail_lines=taxonomy_details,
    )
    logger.add_divider(length=120, border="+", fill="-")

    # Step 3: Download genomes (domain stage)
    genomes_root = dataset_root / genus_name / "genomes"
    already_registered = is_genus_genomes_registered(
        genus_name=genus_name,
        data_root=dataset_root,
    )

    # Initialize state variables
    download_state: str
    registry_state: str
    assemblies_df = None
    registry_details: list[str] = []

    # Case 1: dataset already registered and overwrite disabled
    if already_registered and not overwrite:
        download_state = "skipped"
        registry_state = "dataset_reused"
        registry_details = [
            "  - Existing genomes dataset found",
            "  - Overwrite is disabled",
            "  - No download executed",
        ]

    # Case 2: execute download
    else:
        try:
            assemblies_df, genomes_root = download_genomes_from_genus(
                genus_name=genus_name,
                species_table=species_table,
                data_root=dataset_root,
                assembly_level=assembly_level,
                refseq_categories=refseq_categories,
                source=source,
                parallel=parallel,
                batch_size=batch_size,
                overwrite=overwrite,
            )
            download_state = "executed"

        except ValueError as exc:
            logger.log_with_borders(
                level=logging.INFO, message="[ Download genomes ]", length=120
            )
            raise PipelineError(
                pipeline="genomes.genus",
                stage="download",
                reason=str(exc),
                action="Check taxonomy resolution and species taxids.",
                context={"genus": genus_name},
            ) from exc

        except RuntimeError as exc:
            logger.log_with_borders(
                level=logging.INFO, message="[ Download genomes ]", length=120
            )
            raise PipelineError(
                pipeline="genomes.genus",
                stage="download",
                reason="Genome download execution failed.",
                action="Check NCBI availability or network connectivity.",
                context={"genus": genus_name},
            ) from exc

        # Log download summary
        if assemblies_df.empty:
            registry_state = "registry_not_applicable"
            registry_details = [
                "  - Genome download executed successfully",
                "  - No assemblies matched the retrieval criteria",
                "  - Dataset registration skipped",
            ]
        else:
            try:
                dataset_yaml, assemblies_csv = register_genus_genomes_dataset(
                    genus_name=genus_name,
                    data_root=dataset_root,
                    assemblies=assemblies_df,
                    parameters={
                        "assembly_level": assembly_level,
                        "refseq_categories": refseq_categories,
                        "source": source,
                        "parallel": parallel,
                        "batch_size": batch_size,
                    },
                    overwrite=overwrite,
                )
                registry_state = "dataset_registered"
                registry_details = [
                    "  - Genome download executed successfully",
                    "  - Registry artifacts created",
                    "  - Artifacts :",
                    f"      • Dataset manifest : '{dataset_yaml}'",
                    f"      • Assemblies table : '{assemblies_csv}'",
                ]

            except FileExistsError as exc:
                registry_state = "dataset_reused"
                registry_details = [
                    "  - Genome download executed successfully",
                    "  - Registry artifacts already exist",
                    "  - Existing dataset reused (overwrite disabled)",
                    f"  - Reason: {str(exc)}",
                ]

    assemblies_repr = "N/A" if assemblies_df is None else str(len(assemblies_df))
    render_stage_summary(
        logger=logger,
        title="[ Download genomes ]",
        summary_lines=[
            f"▸ Download state  : '{download_state}'",
            f"▸ Registry state  : '{registry_state}'",
            f"▸ Assemblies      : {assemblies_repr}",
            f"▸ Dataset path    : '{genomes_root}'",
        ],
        detail_lines=registry_details,
    )

    return genomes_root
