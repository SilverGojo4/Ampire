# pylint: disable=
"""
Registry utilities for persisting and accessing
tblastn discovery search artifacts in the Ampire data hierarchy.

This module defines the registry-layer contract for discovery searches.
It is responsible for persisting structured discovery summaries and
scientific fingerprints into stable, reproducible YAML artifacts.
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import Dict

# App Imports
from app.utils.config_utils import load_yaml, write_yaml


# Discovery Registry
def register_discovery(
    *,
    search_id: str,
    data_root: Path,
    discovery_summary: Dict,
    overwrite: bool = False,
) -> Path:
    """
    Register a tblastn discovery dataset artifact.

    This function persists a structured discovery search summary
    into the Ampire data hierarchy as a YAML registry artifact.

    Parameters
    ----------
    search_id : str
        Deterministic search identifier produced from
        `compute_search_id`.
    data_root : pathlib.Path
        Root directory for discovery datasets.
    discovery_summary : dict
        Structured discovery metadata including fingerprint,
        parameters, and produced artifacts.
    overwrite : bool, optional
        Whether to overwrite an existing registry artifact.

    Returns
    -------
    pathlib.Path
        Path to the registered discovery YAML artifact.

    Raises
    ------
    FileExistsError
        If the registry artifact already exists and overwrite
        is False.
    """
    data_root = data_root.resolve()

    # Construct registry paths
    discovery_dir = data_root / search_id
    discovery_yaml_path = discovery_dir / "discovery.yaml"

    # Check existing registration
    if discovery_yaml_path.exists() and not overwrite:
        raise FileExistsError(
            "Discovery registry already exists and overwrite is disabled."
        )

    # Create directory only when registering
    discovery_dir.mkdir(parents=True, exist_ok=True)

    # Persist registry artifact
    write_yaml(
        discovery_yaml_path,
        discovery_summary,
        overwrite=overwrite,
    )

    return discovery_yaml_path


def is_discovery_registered(
    *,
    search_id: str,
    data_root: Path,
) -> bool:
    """
    Determine whether a discovery dataset is registered.

    Parameters
    ----------
    search_id : str
        Deterministic discovery search identifier.
    data_root : pathlib.Path
        Root discovery data directory.

    Returns
    -------
    bool
        True if the registry artifact exists; False otherwise.
    """
    discovery_dir = data_root / search_id
    discovery_yaml = discovery_dir / "discovery.yaml"
    return discovery_yaml.is_file()


def load_discovery(
    *,
    search_id: str,
    data_root: Path,
) -> Dict:
    """
    Load a registered discovery dataset artifact.

    Parameters
    ----------
    search_id : str
        Deterministic discovery search identifier.
    data_root : pathlib.Path
        Root discovery data directory.

    Returns
    -------
    dict
        Discovery registry metadata loaded from YAML.

    Raises
    ------
    FileNotFoundError
        If the discovery dataset has not been registered.
    """
    discovery_dir = data_root / search_id
    discovery_yaml = discovery_dir / "discovery.yaml"

    if not discovery_yaml.is_file():
        raise FileNotFoundError(f"Discovery dataset not registered: '{search_id}'")

    return load_yaml(discovery_yaml)
