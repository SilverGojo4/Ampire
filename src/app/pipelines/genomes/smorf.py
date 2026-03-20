# pylint: disable=too-many-locals, broad-exception-caught, too-many-arguments, too-many-branches, too-many-statements
"""
Genus-level smORF dataset pipeline.

This module orchestrates the end-to-end workflow for building a
reproducible smORF dataset for a given genus, including:

1. Validating the registered genomes dataset
2. Executing smORFinder on each genome FASTA
3. Assembling normalized raw smORF predictions
4. Organizing candidate smORF predictions (legacy smorf outputs)
5. Registering dataset artifacts and provenance
"""

from __future__ import annotations

# Standard Library Imports
import csv
import logging
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional

# App Imports
from app.exceptions import PipelineError
from app.inputs.genus import load_genus_list
from app.processing.genomes.smorf.parser import assemble_smorf_predictions
from app.processing.genomes.smorf.postprocess import (
    convert_predictions_tsv_to_csv,
    rename_smorf_outputs,
    save_smorf_table,
)
from app.processing.genomes.smorf.run_smorf import run_smorf_single_genome
from app.registry.genomes import is_genus_genomes_registered, load_genus_genomes_dataset
from app.registry.smorfs import (
    is_genus_smorfs_registered,
    register_genus_smorfs_dataset,
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
def build_genus_smorfs_batch_dataset(
    *,
    ctx: ExecutionContext,
    inputs: Dict,
    parameters: Dict,
    outputs: Dict,
    execution: Dict,
) -> Dict:
    """
    Build genus-level smORF datasets in batch mode.

    This pipeline orchestrates the execution of the single-genus
    smORF dataset pipeline across multiple genera.

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
        Algorithm and tool parameters controlling smORF prediction.
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
                    "smorfs_root": Optional[str],
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
    dataset_root = Path(outputs["dataset_root"])

    # Parse algorithm parameters
    dsn1_cutoff = parameters.get("dsn1_cutoff", 0.9999)
    dsn2_cutoff = parameters.get("dsn2_cutoff", 0.9999)
    phmm_cutoff = parameters.get("phmm_cutoff", 1e-6)
    smorf_conda_env = parameters.get("smorf_conda_env")

    # Load genus list
    genus_list = load_genus_list(genus_list_path)
    total = len(genus_list)

    # Prepare batch logging context
    pipeline_context = {
        "Input Summary": {
            "Genus list path": str(genus_list_path),
            "Total genera": total,
        },
        "Prediction Parameters": {
            "DSN1 cutoff": dsn1_cutoff,
            "DSN2 cutoff": dsn2_cutoff,
            "pHMM cutoff": phmm_cutoff,
        },
        "Runtime": {"smORFinder env": smorf_conda_env or "current environment"},
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
        title="Genus smORF Batch Pipeline",
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
            smorfs_root = _build_genus_smorfs_dataset(
                ctx=ctx,
                genus_name=genus_name,
                dataset_root=dataset_root,
                dsn1_cutoff=dsn1_cutoff,
                dsn2_cutoff=dsn2_cutoff,
                phmm_cutoff=phmm_cutoff,
                smorf_conda_env=smorf_conda_env,
                overwrite=overwrite,
            )

            results.append(
                {
                    "genus": genus_name,
                    "status": "success",
                    "smorfs_root": str(smorfs_root),
                    "stage": None,
                    "reason": None,
                }
            )
            success += 1

        except PipelineError as exc:
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
                    "smorfs_root": None,
                }
            )
            failed += 1

        except Exception as exc:
            fallback_error = PipelineError(
                pipeline="genomes.smorf",
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
                    "smorfs_root": None,
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
        title="Genus smORF Batch Pipeline",
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
def _build_genus_smorfs_dataset(
    *,
    ctx: ExecutionContext,
    genus_name: str,
    dataset_root: Path,
    dsn1_cutoff: float,
    dsn2_cutoff: float,
    phmm_cutoff: float,
    smorf_conda_env: Optional[str] = None,
    overwrite: bool = False,
) -> Optional[Path]:
    """
    Build a genus-level smORF dataset.

    Parameters
    ----------
    ctx : ExecutionContext
        Runtime execution context inherited from the batch pipeline.
    genus_name : str
        Scientific genus name.
    dataset_root : pathlib.Path
        Root directory under which all genus-level dataset artifacts
        will be created.
    dsn1_cutoff : float
        DSN1 probability cutoff.
    dsn2_cutoff : float
        DSN2 probability cutoff.
    phmm_cutoff : float
        Profile HMM significance cutoff.
    overwrite : bool, optional
        Whether to overwrite existing registered smORF datasets.

    Returns
    -------
    Optional[pathlib.Path]
        Path to the registered genus-level smORF dataset directory.

    Raises
    ------
    PipelineError
        If any pipeline stage fails.
    """
    # Initialize logger
    logger = ctx.logger

    # Step 1: Validate genomes dataset dependency
    if not is_genus_genomes_registered(
        genus_name=genus_name,
        data_root=dataset_root,
    ):
        logger.log_with_borders(
            level=logging.INFO,
            message="[ Validate genomes dataset ]",
            length=120,
        )
        raise PipelineError(
            pipeline="genomes.smorf",
            stage="inputs",
            reason="Required genomes dataset is not registered.",
            action="Run the 'genomes.genus' pipeline first.",
            context={"genus": genus_name},
        )

    genomes_metadata, assemblies = load_genus_genomes_dataset(
        genus_name=genus_name,
        data_root=dataset_root,
    )

    fasta_dir = dataset_root / genus_name / "genomes" / "fasta"
    smorfs_root = dataset_root / genus_name / "smorfs"
    raw_root = smorfs_root / "raw"
    pool_root = smorfs_root / "pool"
    confidence_root = smorfs_root / "conf"

    dependency_details = [
        f"  - Genomes dataset : '{dataset_root / genus_name / 'genomes'}'",
        f"  - FASTA directory : '{fasta_dir}'",
        f"  - Assemblies      : {len(assemblies)}",
    ]
    render_stage_summary(
        logger=logger,
        title="[ Validate genomes dataset ]",
        summary_lines=[
            "▸ Dependency state : 'available'",
            f"▸ Dataset path     : '{dataset_root / genus_name / 'genomes'}'",
        ],
        detail_lines=dependency_details,
    )
    logger.add_divider(length=120, border="+", fill="-")

    # Step 2: Check reuse policy for smORF dataset
    already_registered = is_genus_smorfs_registered(
        genus_name=genus_name,
        data_root=dataset_root,
    )

    if already_registered and not overwrite:
        render_stage_summary(
            logger=logger,
            title="[ Register smORF dataset ]",
            summary_lines=[
                "▸ Registry state : 'dataset_reused'",
                f"▸ Dataset path   : '{smorfs_root}'",
            ],
            detail_lines=[
                "  - Existing smORF dataset found",
                "  - Overwrite is disabled",
                "  - No execution performed",
            ],
        )
        return smorfs_root

    # Step 3: Prepare output directories
    if overwrite and smorfs_root.exists():
        shutil.rmtree(smorfs_root)

    raw_root.mkdir(parents=True, exist_ok=True)
    pool_root.mkdir(parents=True, exist_ok=True)
    confidence_root.mkdir(parents=True, exist_ok=True)

    # Step 4: Execute smORFinder genome by genome
    total = 0
    success_count = 0
    failed_genomes: list[str] = []

    pool_aggregate_records: list[dict[str, str]] = []
    conf_aggregate_records: list[dict[str, str]] = []

    for _, row in assemblies.iterrows():
        total += 1

        local_path = Path(row["local_filename"])

        genome_filename = local_path.name
        genome_id = local_path.stem
        fasta_path = fasta_dir / genome_filename

        if not fasta_path.is_file():
            failed_genomes.append(genome_id)
            continue

        raw_genome_dir = raw_root / genome_id
        pool_genome_dir = pool_root / genome_id
        confidence_genome_dir = confidence_root / genome_id

        raw_genome_dir.mkdir(parents=True, exist_ok=True)
        pool_genome_dir.mkdir(parents=True, exist_ok=True)
        confidence_genome_dir.mkdir(parents=True, exist_ok=True)

        # 4.1 Run smORFinder
        ok = run_smorf_single_genome(
            fasta_path=fasta_path,
            output_dir=raw_genome_dir,
            dsn1_cutoff=dsn1_cutoff,
            dsn2_cutoff=dsn2_cutoff,
            phmm_cutoff=phmm_cutoff,
            smorf_conda_env=smorf_conda_env,
            overwrite=overwrite,
        )

        if not ok:
            failed_genomes.append(genome_id)
            continue

        # 4.2 Assemble normalized raw predictions
        records = assemble_smorf_predictions(
            genome_dir=raw_genome_dir,
        )

        # 4.3 Assemble final prediction table
        save_smorf_table(
            records=records,
            output_path=pool_genome_dir / "predictions.tsv",
        )

        for rec in records:
            rec_copy = dict(rec)
            rec_copy["genome_id"] = genome_id
            pool_aggregate_records.append(rec_copy)

        # 4.4 Copy key small-ORF files
        _copy_smorf_pool_artifacts(
            src_dir=raw_genome_dir,
            dst_dir=pool_genome_dir,
        )
        _copy_smorf_conf_artifacts(
            src_dir=raw_genome_dir,
            dst_dir=confidence_genome_dir,
        )

        # 4.5 Collect high-confidence predictions from raw TSV
        conf_tsv_candidates = list(raw_genome_dir.glob("*.tsv"))
        if len(conf_tsv_candidates) != 1:
            raise PipelineError(
                pipeline="genomes.smorf",
                stage="aggregate_conf",
                reason="Expected exactly one raw TSV output for high-confidence predictions.",
                action="Inspect smORFinder raw outputs for this genome.",
                context={
                    "genus": genus_name,
                    "genome_id": genome_id,
                    "raw_genome_dir": str(raw_genome_dir),
                    "tsv_count": len(conf_tsv_candidates),
                },
            )

        conf_tsv_path = conf_tsv_candidates[0]

        with open(conf_tsv_path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row_dict in reader:
                row_dict["genome_id"] = genome_id
                conf_aggregate_records.append(row_dict)

        # 4.6 Clean + normalize dataset outputs
        rename_smorf_outputs(
            genome_output_dir=pool_genome_dir,
        )
        convert_predictions_tsv_to_csv(
            genome_output_dir=pool_genome_dir,
        )
        rename_smorf_outputs(
            genome_output_dir=confidence_genome_dir,
        )
        convert_predictions_tsv_to_csv(
            genome_output_dir=confidence_genome_dir,
        )

        success_count += 1

    # Step 5: Save aggregate prediction tables
    pool_aggregate_csv = pool_root / "ALL_predictions.csv"
    conf_aggregate_csv = confidence_root / "ALL_predictions.csv"

    if pool_aggregate_records:
        pool_fieldnames = ["genome_id"] + [
            k for k in pool_aggregate_records[0].keys() if k != "genome_id"
        ]
        with open(pool_aggregate_csv, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=pool_fieldnames)
            writer.writeheader()
            writer.writerows(pool_aggregate_records)
    else:
        logger.warning(f"No pool aggregate records collected for genus '{genus_name}'.")

    if conf_aggregate_records:
        conf_fieldnames = ["genome_id"] + [
            k for k in conf_aggregate_records[0].keys() if k != "genome_id"
        ]
        with open(conf_aggregate_csv, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=conf_fieldnames)
            writer.writeheader()
            writer.writerows(conf_aggregate_records)
    else:
        logger.warning(f"No conf aggregate records collected for genus '{genus_name}'.")

    render_stage_summary(
        logger=logger,
        title="[ Execute smORFinder ]",
        summary_lines=[
            f"▸ Total genomes              : {total}",
            f"▸ Successfully done          : {success_count}",
            f"▸ Failed genomes             : {len(failed_genomes)}",
            f"▸ Pool aggregate records     : {len(pool_aggregate_records)}",
            f"▸ Conf aggregate records     : {len(conf_aggregate_records)}",
            f"▸ Raw output path            : '{raw_root}'",
            f"▸ Candidate pool output path : '{pool_root}'",
            f"▸ High confidence path       : '{confidence_root}'",
        ],
        detail_lines=(
            [f"  - Failed list : {', '.join(failed_genomes)}"]
            if failed_genomes
            else ["  - All genomes processed successfully"]
        ),
    )
    logger.add_divider(length=120, border="+", fill="-")

    # Step 6: Register dataset
    try:
        dataset_yaml = register_genus_smorfs_dataset(
            genus_name=genus_name,
            data_root=dataset_root,
            parameters={
                "dsn1_cutoff": dsn1_cutoff,
                "dsn2_cutoff": dsn2_cutoff,
                "phmm_cutoff": phmm_cutoff,
                "genomes_dataset": genomes_metadata.get("dataset", {}).get("name"),
            },
            overwrite=overwrite,
        )
        registry_state = "dataset_registered"
        registry_details = [
            "  - smORF prediction executed successfully",
            "  - Registry artifacts created",
            "  - Artifacts :",
            f"      • Dataset manifest : '{dataset_yaml}'",
            f"      • Raw dataset      : '{raw_root}'",
            f"      • Candidate pool   : '{pool_root}'",
            f"      • High confidence  : '{confidence_root}'",
        ]
    except FileExistsError as exc:
        registry_state = "dataset_reused"
        registry_details = [
            "  - Registry artifacts already exist",
            "  - Existing dataset reused (overwrite disabled)",
            f"  - Reason: {str(exc)}",
        ]

    render_stage_summary(
        logger=logger,
        title="[ Register smORF dataset ]",
        summary_lines=[
            f"▸ Registry state : '{registry_state}'",
            f"▸ Dataset path   : '{smorfs_root}'",
        ],
        detail_lines=registry_details,
    )

    return smorfs_root


def _copy_smorf_pool_artifacts(
    *,
    src_dir: Path,
    dst_dir: Path,
) -> None:
    """
    Copy smORF sequence artifacts between directories.

    Parameters
    ----------
    src_dir : pathlib.Path
        Source smORF output directory.
    dst_dir : pathlib.Path
        Destination directory.
    """

    artifact_files = [
        "prodigal.small.faa",
        "prodigal.small.ffn",
        "prodigal.small.gff",
    ]

    for fname in artifact_files:

        src = src_dir / "tmp" / fname
        dst = dst_dir / fname

        if src.exists() and not dst.exists():
            shutil.copy(src, dst)


def _copy_smorf_conf_artifacts(
    *,
    src_dir: Path,
    dst_dir: Path,
) -> None:
    """
    Copy final smORF outputs from raw directory to conf directory.

    Parameters
    ----------
    src_dir : pathlib.Path
        Source smORF output directory.
    dst_dir : pathlib.Path
        Destination directory.
    """
    suffixes = {".faa", ".ffn", ".gff", ".tsv"}

    for file in src_dir.iterdir():
        if file.is_file() and file.suffix in suffixes:
            shutil.copy(file, dst_dir / file.name)
