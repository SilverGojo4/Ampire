# pylint: disable=
"""
Execution context for Ampire pipeline runs.

This module defines a lightweight runtime context object
that is passed through the pipeline execution stack.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass


@dataclass(frozen=True)
class ExecutionContext:
    """
    Runtime execution context for a single pipeline run.

    Attributes
    ----------
    pipeline_name : str
        Identifier of the pipeline being executed.
    run_id : str
        Unique execution identifier for this run.
    logger : logging.Logger
        Logger instance bound to this execution.
    """

    pipeline_name: str
    run_id: str
    logger: logging.Logger | logging.LoggerAdapter
