# pylint: disable=missing-function-docstring
"""
Ampire command-line interface.

This module is responsible for:
- Parsing CLI arguments
- Loading pipeline configuration (YAML)
- Dispatching execution to the correct pipeline
"""

from __future__ import annotations

# Standard Library
import argparse
import importlib
import logging
from pathlib import Path
from typing import Callable, Dict

# App imports
from app.exceptions import PipelineError
from app.runtime.context import ExecutionContext
from app.runtime.logging.initializer import initialize_pipeline_logger
from app.utils.config_utils import load_yaml

# Pipeline registry (dispatcher table)
# Format: "pipeline_name": "module.path:function_name"
PIPELINES = {
    "genomes.genus": "app.pipelines.genomes.genus:build_genus_genomes_batch_dataset",
    "genomes.blastdb_global": "app.pipelines.genomes.blastdb_global:build_global_blastdb_pipeline",
    # future:
    # "genomes.species": "app.pipelines.genomes.species:build_species_genomes_dataset",
    # "proteins.genus": "app.pipelines.proteins.genus:build_genus_proteins_dataset",
}


# Public CLI entrypoint
def main() -> None:
    """
    Entry point for the Ampire command-line interface.

    This function parses CLI arguments and dispatches execution
    to the appropriate pipeline runner.
    """
    parser = argparse.ArgumentParser(
        prog="ampire",
        description="Ampire pipeline runner",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # ampire run
    run_parser = subparsers.add_parser(
        "run",
        help="Run a pipeline",
    )
    run_parser.add_argument(
        "--pipeline",
        required=True,
        choices=PIPELINES.keys(),
        help="Pipeline identifier (e.g. genomes.genus)",
    )
    run_parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to pipeline YAML config",
    )

    args = parser.parse_args()

    if args.command == "run":
        _run_pipeline(
            pipeline_name=args.pipeline,
            config_path=args.config,
        )


# Pipeline execution logic
def _run_pipeline(
    *,
    pipeline_name: str,
    config_path: Path,
) -> None:
    """
    Execute a registered pipeline from the command-line interface.

    This internal helper function serves as the boundary between
    the Ampire CLI and pipeline execution logic. It loads pipeline
    configuration, resolves the selected pipeline entry point,
    and invokes the pipeline with standardized execution arguments.

    Pipeline errors are intercepted and rendered as user-friendly
    CLI messages rather than propagating as uncaught exceptions.

    Parameters
    ----------
    pipeline_name : str
        Identifier of the pipeline to execute (e.g., "genomes.genus").
    config_path : pathlib.Path
        Path to the pipeline YAML configuration file.

    Returns
    -------
    None
    """
    # Lazy pipeline loading
    pipeline_path = PIPELINES[pipeline_name]
    pipeline_fn = _load_pipeline_callable(pipeline_path)

    # Load config
    config = _parse_pipeline_config(config_path)

    # inputs = _normalize_pipeline_inputs(config["inputs"])
    inputs = config["inputs"]
    parameters = config["parameters"]
    outputs = config["outputs"]
    execution = config["execution"]

    ctx = _create_execution_context(pipeline_name=pipeline_name)

    # Execute pipeline
    try:
        pipeline_fn(
            ctx=ctx,
            inputs=inputs,
            parameters=parameters,
            outputs=outputs,
            execution=execution,
        )
    except PipelineError as exc:
        _print_pipeline_error(exc)


# Internal Helpers
def _load_pipeline_callable(pipeline_path: str) -> Callable:
    """
    Resolve and load a pipeline entry-point callable.

    Parameters
    ----------
    pipeline_path : str
        Fully qualified pipeline reference in the format
        "module.path:function_name".

    Returns
    -------
    Callable
        Resolved pipeline entry-point function.

    Raises
    ------
    PipelineError
        If the pipeline path format is invalid, the module
        cannot be imported, or the function is not found.
    """
    # Parse pipeline reference
    try:
        module_path, func_name = pipeline_path.split(":")
    except ValueError as exc:
        raise PipelineError(
            pipeline="execution",
            stage="resolve",
            reason="Invalid pipeline path format",
            action=(
                "Pipeline path must follow the format " "'module.path:function_name'."
            ),
            context={"pipeline_path": pipeline_path},
        ) from exc

    # Import pipeline module
    module = importlib.import_module(module_path)

    # Resolve entry-point callable
    try:
        return getattr(module, func_name)
    except AttributeError as exc:
        raise PipelineError(
            pipeline="execution",
            stage="resolve",
            reason="Pipeline entry-point function not found",
            action="Ensure the function name exists in the target module.",
            context={
                "module": module_path,
                "function": func_name,
            },
        ) from exc


def _parse_pipeline_config(config_path: Path) -> Dict[str, Dict]:
    """
    Load and split a pipeline configuration file into logical sections.

    Parameters
    ----------
    config_path : pathlib.Path
        Path to the pipeline YAML configuration file.

    Returns
    -------
    Dict[str, Dict]
        Dictionary containing the following sections:
        - inputs
        - parameters
        - outputs
        - execution
    """
    config = load_yaml(config_path)

    return {
        "inputs": config.get("inputs", {}),
        "parameters": config.get("parameters", {}),
        "outputs": config.get("outputs", {}),
        "execution": config.get("execution", {}),
    }


def _create_execution_context(*, pipeline_name: str) -> ExecutionContext:
    """
    Create execution context for a single pipeline run.

    Parameters
    ----------
    pipeline_name : str
        Logical pipeline identifier.

    Returns
    -------
    ExecutionContext
        Initialized execution context instance.
    """
    # Initialize logging backend
    logger = initialize_pipeline_logger(
        pipeline_name=pipeline_name,
        logging_config_path=Path("configs/logging/general.json"),
        logs_root=Path("logs"),
    )

    logger.info("[ 'Pipeline Initialization Summary' ]")
    logger.log_pipeline_initialization(
        project_name="Ampire",
        line_width=120,
    )
    logger.add_spacer(level=logging.INFO, lines=1)

    # Bind execution context
    return ExecutionContext(
        pipeline_name=pipeline_name,
        run_id=logger.name if hasattr(logger, "name") else "unknown",
        logger=logger,
    )


def _print_pipeline_error(exc: PipelineError) -> None:
    """
    Render pipeline errors in a human-readable CLI format.
    """
    print()
    print("[PIPELINE ERROR]")
    print(f"Pipeline: {exc.pipeline}")
    print(f"Stage   : {exc.stage}")

    if exc.context and "genus" in exc.context:
        print(f"Genus   : {exc.context['genus']}")

    print()
    print("Reason:")
    print(f"  {exc.reason}")

    if exc.action:
        print()
        print("Action:")
        print(f"  {exc.action}")

    print()


if __name__ == "__main__":
    main()
