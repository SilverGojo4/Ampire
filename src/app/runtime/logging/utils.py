# pylint: disable=
"""
Utility helpers for pipeline runtime logging.

This module provides small, reusable formatting helpers
used by logging renderers and pipeline execution code
to ensure consistent, human-readable output.
"""

from __future__ import annotations

# Standard Library Imports
from datetime import datetime, timedelta, timezone


def format_duration(duration: timedelta) -> str:
    """
    Format a timedelta into a concise, human-readable
    duration string.

    Parameters
    ----------
    duration : datetime.timedelta
        Duration to be formatted.

    Returns
    -------
    str
        Human-readable duration string.
    """
    total_seconds = int(duration.total_seconds())

    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    parts = []

    if hours:
        parts.append(f"{hours}h")
    if minutes:
        parts.append(f"{minutes}m")
    if seconds or not parts:
        parts.append(f"{seconds}s")

    return " ".join(parts)


def format_timestamp(dt: datetime) -> str:
    """
    Format a timezone-aware datetime into a log-friendly
    timestamp string.

    Parameters
    ----------
    dt : datetime.datetime
        Timezone-aware datetime (UTC recommended).

    Returns
    -------
    str
        Formatted timestamp, e.g. "2025-12-02 12:27:25"
    """
    if dt.tzinfo is None:
        raise ValueError("format_timestamp requires a timezone-aware datetime")

    return dt.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
