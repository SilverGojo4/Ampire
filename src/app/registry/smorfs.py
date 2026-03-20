# pylint: disable=
"""
Registry utilities for persisting and accessing
genus-level smORF dataset artifacts in the Ampire data hierarchy.
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import Dict

# Third-Party Library Imports
import pandas as pd

# App Imports
from app.utils.config_utils import load_yaml, write_yaml


# smORF Registry
def register_genus_smorfs_dataset(
    *,
    genus_name: str,
    data_root: Path,
    parameters: Dict,
    overwrite: bool = False,
) -> Path:
    """
    Register smORF dataset artifacts for a given genus.

    This function persists the smORF discovery results
    into the Ampire dataset hierarchy by writing a
    dataset manifest describing the dataset and its
    associated artifacts.

    Parameters
    ----------
    genus_name : str
        Scientific genus name.
    data_root : pathlib.Path
        Root directory where genus datasets are stored.
    parameters : dict
        Parameters used during smORF prediction.
    overwrite : bool, optional
        Whether to overwrite an existing registered dataset.

    Returns
    -------
    pathlib.Path
        Path to the dataset manifest file.

    Raises
    ------
    FileExistsError
        If dataset artifacts already exist and overwrite is False.
    """

    genus_root = data_root / genus_name
    smorfs_dir = genus_root / "smorfs"

    dataset_yaml = smorfs_dir / "dataset.yaml"

    if dataset_yaml.exists() and not overwrite:
        raise FileExistsError(
            "smORF dataset artifacts already exist and overwrite is disabled."
        )

    smorfs_dir.mkdir(parents=True, exist_ok=True)

    dataset_manifest = {
        "dataset": {
            "name": "smorfs",
            "scope": "genus",
            "genus": genus_name,
        },
        "source": {
            "tool": "smORFinder",
        },
        "parameters": parameters,
        "provenance": {
            "registered_at": pd.Timestamp.utcnow().isoformat(),
            "registered_by": "register_genus_smorfs_dataset",
        },
        "artifacts": {
            "raw_dir": "raw/",
            "pool_dir": "pool/",
            "conf_dir": "conf/",
        },
    }

    write_yaml(
        dataset_yaml,
        dataset_manifest,
        overwrite=overwrite,
    )

    return dataset_yaml


def is_genus_smorfs_registered(
    genus_name: str,
    data_root: Path,
) -> bool:
    """
    Determine whether a smORF dataset has been registered.

    Parameters
    ----------
    genus_name : str
        Scientific genus name.
    data_root : pathlib.Path
        Root directory for dataset storage.

    Returns
    -------
    bool
        True if the dataset manifest exists.
    """

    dataset_yaml = data_root / genus_name / "smorfs" / "dataset.yaml"

    return dataset_yaml.is_file()


def load_genus_smorfs_dataset(
    genus_name: str,
    data_root: Path,
) -> Dict:
    """
    Load registered smORF dataset metadata.

    Parameters
    ----------
    genus_name : str
        Scientific genus name.
    data_root : pathlib.Path
        Root directory where datasets are stored.

    Returns
    -------
    dict
        Dataset manifest contents.

    Raises
    ------
    FileNotFoundError
        If the smORF dataset has not been registered.
    """

    dataset_yaml = data_root / genus_name / "smorfs" / "dataset.yaml"

    if not dataset_yaml.exists():
        raise FileNotFoundError(
            f"smORF dataset not registered for genus '{genus_name}'."
        )

    return load_yaml(dataset_yaml)
