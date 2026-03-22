# pylint: disable=too-few-public-methods
"""
smORFinder error definitions.

This module defines structured exception classes used during
smORF prediction execution. These exceptions enable precise
classification of failure modes and allow upstream pipelines
to handle errors in a fault-tolerant manner.
"""

from __future__ import annotations


# Base class
class SmorfError(Exception):
    """Base class for all smORFinder-related errors."""


# Execution-level errors
class SmorfExecutionError(SmorfError):
    """
    Generic smORFinder execution failure.

    Raised when the external `smorf` command fails with
    a non-zero return code and cannot be further classified.
    """

    def __init__(self, message: str, *, cmd=None, stderr=None):
        super().__init__(message)
        self.cmd = cmd
        self.stderr = stderr


class SmorfRuntimeCrashError(SmorfExecutionError):
    """
    Low-level runtime crash (e.g., segmentation fault, invalid pointer).

    Typically originates from underlying tools such as Prodigal.
    """


class SmorfInternalError(SmorfExecutionError):
    """
    Internal Python exception within smORFinder.

    Indicates a traceback occurred inside the smORFinder pipeline.
    """


class SmorfDependencyError(SmorfExecutionError):
    """
    Error related to missing or misconfigured dependencies.

    Example:
    - `smorf` executable not found
    - Conda environment issues
    """


# GFF / Prodigal errors
class EmptyGFFError(SmorfError):
    """
    No genes predicted (empty GFF output).
    """


class MalformedGFFError(SmorfError):
    """
    GFF file exists but has invalid format.
    """


class GFFParsingError(SmorfError):
    """
    Error occurred during parsing of GFF file.

    Example:
    - IndexError in smorfinder.prodigal
    """


# Output validation errors
class MissingOutputError(SmorfError):
    """
    Expected output files are missing.
    """


class IncompleteOutputError(SmorfError):
    """
    Output exists but is incomplete or corrupted.
    """


# Input validation errors
class InvalidFASTAError(SmorfError):
    """
    Input FASTA file is invalid or malformed.
    """
