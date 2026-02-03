# pylint: disable=too-many-arguments
"""
Structured log renderers for pipeline runtime output.

This module provides reusable helpers for rendering
well-formatted, sectioned log output from structured
logging contexts.
"""

from __future__ import annotations

# Standard Library Imports
import logging
from typing import Dict

# App Imports
from app.exceptions import PipelineError


def render_batch_header(
    *,
    logger,
    title: str,
    context: Dict[str, Dict[str, object]],
) -> None:
    """
    Render a batch-level pipeline header from a structured
    logging context.

    Parameters
    ----------
    logger
        Pipeline logger instance.
    title : str
        Human-readable pipeline title.
    context : dict
        Nested mapping of section titles to key-value pairs
        describing batch-level configuration and inputs.

        Example :
        {
            "Input Summary": {
                "Total genera": 128,
                "Dataset root": "data/interim",
            },
            "Retrieval Parameters": {
                "Source": "refseq",
                "Assembly level": "complete,chromosome",
            },
        }

    Returns
    -------
    None
    """
    logger.info(f"[ '{title}' ]")
    logger.add_divider(length=120, border="+", fill="=")

    first = True
    for section, items in context.items():

        # Separate sections visually (but lightly)
        if not first:
            logger.add_divider(length=120, border="║", fill="-")
        first = False

        lines = [f"[ {section} ]"]
        for key, value in items.items():
            lines.append(f"▸ {key:<20} : {value}")

        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(lines),
            border="║",
            length=120,
        )

    logger.add_divider(length=120, border="+", fill="=")
    logger.add_spacer(level=logging.INFO, lines=1)


def render_batch_footer(
    *,
    logger,
    title: str,
    context: Dict[str, Dict[str, object]],
) -> None:
    """
    Render a batch-level pipeline execution footer from a
    structured result context.

    Parameters
    ----------
    logger
        Pipeline logger instance.
    title : str
        Human-readable pipeline title.
    context : dict
        Nested mapping of section titles to key-value pairs
        describing execution results and timing.

        Example :
        {
            "Results": {
                "Total": 128,
                "Success": 120,
                "Failed": 6,
                "Skipped": 2,
            },
            "Timing": {
                "Started at": "2026-02-02T10:31:00",
                "Finished at": "2026-02-02T10:42:18",
                "Duration": "11m 18s",
            },
        }

    Returns
    -------
    None
    """
    logger.info(f"[ '{title}' — Summary ]")
    logger.add_divider(length=120, border="+", fill="=")

    first = True
    for section, items in context.items():

        if not first:
            logger.add_divider(length=120, border="║", fill="-")
        first = False

        lines = [f"[ {section} ]"]
        for key, value in items.items():
            lines.append(f"▸ {key:<15} : {value}")

        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(lines),
            border="║",
            length=120,
        )

    logger.add_divider(length=120, border="+", fill="=")
    logger.add_spacer(level=logging.INFO, lines=1)


def render_pipeline_error(
    *,
    logger,
    error: PipelineError,
) -> None:
    """
    Render a pipeline execution failure from a PipelineError.

    Parameters
    ----------
    logger
        Pipeline logger instance.
    error : PipelineError
        Pipeline error carrying structured failure information.

    Returns
    -------
    None
    """
    logger.log_with_borders(
        level=logging.ERROR,
        message=(
            "▸ Status : 'Failed'\n"
            f"▸ Stage  : '{error.stage}'\n"
            f"▸ Reason : {error.reason}\n"
            f"▸ Action : {error.action}"
        ),
        length=120,
    )


def render_stage_summary(
    *,
    logger,
    title: str | None = None,
    summary_lines: list[str],
    detail_lines: list[str] | None = None,
    level: int = logging.INFO,
    length: int = 120,
) -> None:
    """
    Render a structured summary for a single pipeline stage.

    Parameters
    ----------
    logger
        Pipeline logger instance.
    title : str, optional
        Stage title (e.g. "[ Download genomes ]").
    summary_lines : list[str]
        High-level summary lines (states, counts, paths).
    detail_lines : list[str], optional
        Detailed bullet lines (artifacts, reasons, actions).
    level : int, optional
        Logging level.
    length : int, optional
        Border width.

    Returns
    -------
    None
    """
    lines: list[str] = []

    if title:
        lines.append(title)

    lines.extend(summary_lines)

    if detail_lines:
        lines.append("▸ Details")
        lines.extend(detail_lines)

    logger.log_with_borders(
        level=level,
        message="\n".join(lines),
        length=length,
    )
