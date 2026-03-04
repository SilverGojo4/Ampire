# pylint: disable=
"""
Search fingerprint utilities for the tblastn discovery stage.

This module defines the scientific identity of a tblastn search and
provides deterministic fingerprinting utilities used for registry
tracking, reproducibility guarantees, and cache-safe execution within
the Ampire pipeline.
"""

from __future__ import annotations

# Standard Library Imports
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

# Parameters that define scientific search identity
SCIENTIFIC_PARAMETER_KEYS: List[str] = [
    "evalue",
    "word_size",
    "seg",
    "soft_masking",
    "comp_based_stats",
    "use_sw_tback",
    "max_target_seqs",
    "max_hsps",
    "qcov_hsp_perc",
]


# Public API
def compute_search_fingerprint(
    *,
    query_fasta: Path,
    db_prefix: Path,
    outfmt_fields: List[str],
    parameters: Dict[str, Any],
    blast_version: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Construct a fingerprint describing the scientific identity of a search.

    The fingerprint captures all inputs and parameters that influence
    biological or computational search results, excluding runtime-only
    parameters (e.g., thread count).

    Parameters
    ----------
    query_fasta : pathlib.Path
        Query protein FASTA used for tblastn search.
    db_prefix : pathlib.Path
        BLAST database prefix path.
    outfmt_fields : list[str]
        Ordered BLAST outfmt 6 fields defining output schema.
    parameters : dict[str, Any]
        tblastn execution parameters.
    blast_version : str, optional
        BLAST+ version string used for execution.

    Returns
    -------
    dict[str, Any]
        Fingerprint dictionary representing the deterministic
        scientific configuration of the search.
    """

    query_fasta = Path(query_fasta)
    db_prefix = Path(db_prefix)

    fingerprint: Dict[str, Any] = {
        "tool": "tblastn",
        "task": parameters.get("task", "tblastn"),
        "inputs": {
            "query_fasta": {
                "path": str(query_fasta.resolve()),
                "sha256": _file_sha256(query_fasta),
            },
            "db_prefix": {
                "path": str(db_prefix.resolve()),
                "signature": _blastdb_signature(db_prefix),
            },
        },
        "outfmt": {
            "type": 6,
            "fields": outfmt_fields,
        },
        "parameters": _extract_scientific_parameters(parameters),
    }

    if blast_version:
        fingerprint["environment"] = {"blast_version": blast_version}

    return fingerprint


def compute_search_id(fingerprint: Dict[str, Any]) -> str:
    """
    Compute a deterministic search identifier from a fingerprint.

    Parameters
    ----------
    fingerprint : dict[str, Any]
        Search fingerprint dictionary.

    Returns
    -------
    str
        Short deterministic search identifier derived from SHA256.
    """
    payload = _stable_json_dumps(fingerprint).encode("utf-8")

    # Shortened ID (16 hex chars) remains collision-safe for registry usage
    return hashlib.sha256(payload).hexdigest()[:16]


# Internal Helpers
def _stable_json_dumps(obj: Any) -> str:
    """
    Serialize an object into a canonical JSON string for hashing.

    The serialization guarantees:
    - sorted keys
    - no whitespace
    - UTF-8 safe encoding

    Parameters
    ----------
    obj : Any
        JSON-serializable object.

    Returns
    -------
    str
        Deterministic JSON string representation.
    """
    return json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )


def _file_sha256(path: Path) -> str:
    """
    Compute SHA256 hash of file content.

    Parameters
    ----------
    path : pathlib.Path
        Path to the file whose contents are hashed.

    Returns
    -------
    str
        Hexadecimal SHA256 digest.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    """
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(f"File not found: '{path}'")

    hasher = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            hasher.update(chunk)

    return hasher.hexdigest()


def _blastdb_signature(db_prefix: Path) -> Dict[str, Any]:
    """
    Construct a lightweight filesystem signature for a BLAST database.

    Parameters
    ----------
    db_prefix : pathlib.Path
        BLAST database prefix path. All files matching
        "<prefix>.*" within the parent directory are included
        in the signature.

    Returns
    -------
    dict[str, Any]
        Dictionary describing database identity based on
        filename, file size, and modification time.

    Raises
    ------
    FileNotFoundError
        If the database directory does not exist.
    """
    db_prefix = Path(db_prefix)

    files = sorted(db_prefix.parent.glob(db_prefix.name + ".*"))

    return {
        "files": [
            {
                "name": f.name,
                "size": f.stat().st_size,
                "mtime": int(f.stat().st_mtime),
            }
            for f in files
            if f.is_file()
        ]
    }


def _extract_scientific_parameters(
    parameters: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Extract parameters that influence scientific search results.

    Parameters
    ----------
    parameters : dict[str, Any]
        Full parameter dictionary used for tblastn execution.

    Returns
    -------
    dict[str, Any]
        Subset of parameters contributing to the scientific identity
        of the search.
    """
    scientific_params: Dict[str, Any] = {}

    for key in SCIENTIFIC_PARAMETER_KEYS:
        if key in parameters:
            scientific_params[key] = parameters[key]

    return scientific_params
