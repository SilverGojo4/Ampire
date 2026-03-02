# pylint: disable=
"""
Registry utilities for persisting and accessing
BLAST database artifacts in the Ampire data hierarchy.

This module defines the registry-layer contract for BLAST databases.
It is responsible for persisting structured BLAST DB build summaries
produced by the processing layer into stable, reproducible YAML
artifacts.
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import Dict

# App Imports
from app.utils.config_utils import load_yaml, write_yaml


# BLAST DB Registry
def register_blastdb(
    *,
    db_name: str,
    data_root: Path,
    blastdb_summary: Dict,
    overwrite: bool = False,
) -> Path:
    """
    Register a BLAST database artifact.

    This function persists a structured BLAST database build summary
    into the Ampire data hierarchy as a YAML registry artifact.

    Parameters
    ----------
    db_name : str
        Logical name of the BLAST database (e.g., "ampire_global").
    data_root : pathlib.Path
        Root directory for data storage (e.g., data/processed).
    blastdb_summary : dict
        BLAST database build summary returned by
        `build_global_blast_db`.
    overwrite : bool, optional
        Whether to overwrite an existing BLAST DB registry artifact.
        If False and the registry already exists, an exception is raised.

    Returns
    -------
    pathlib.Path
        Path to the registered BLAST database YAML artifact.

    Raises
    ------
    FileExistsError
        If the BLAST DB registry artifact already exists and
        overwrite is False.
    """
    data_root = data_root.resolve()

    # Construct registry paths
    blastdb_dir = data_root / db_name
    blastdb_yaml_path = blastdb_dir / "blastdb.yaml"

    # Check existing registration
    if blastdb_yaml_path.exists() and not overwrite:
        raise FileExistsError(
            "BLAST DB registry already exists and overwrite is disabled."
        )

    # Create directories (only now we touch filesystem)
    blastdb_dir.mkdir(parents=True, exist_ok=True)

    # Persist registry artifact
    write_yaml(
        blastdb_yaml_path,
        blastdb_summary,
        overwrite=overwrite,
    )

    return blastdb_yaml_path


def is_blastdb_registered(
    *,
    db_name: str,
    data_root: Path,
) -> bool:
    """
    Determine whether a BLAST database is registered.

    Parameters
    ----------
    db_name : str
        Logical name of the BLAST database.
    data_root : pathlib.Path
        Root directory for data storage.

    Returns
    -------
    bool
        True if the BLAST DB registry artifact exists; False otherwise.
    """
    blastdb_dir = data_root / db_name
    blastdb_yaml = blastdb_dir / "blastdb.yaml"
    return blastdb_yaml.is_file()


def load_blastdb(
    *,
    db_name: str,
    data_root: Path,
) -> Dict:
    """
    Load a registered BLAST database artifact.

    Parameters
    ----------
    db_name : str
        Logical name of the BLAST database.
    data_root : pathlib.Path
        Root directory for data storage.

    Returns
    -------
    dict
        BLAST database registry metadata loaded from YAML.

    Raises
    ------
    FileNotFoundError
        If the BLAST database has not been registered.
    """
    blastdb_dir = data_root / db_name
    blastdb_yaml = blastdb_dir / "blastdb.yaml"

    if not blastdb_yaml.is_file():
        raise FileNotFoundError(f"BLAST DB not registered: '{db_name}'")

    return load_yaml(blastdb_yaml)
