# pylint: disable=
"""
Registry utilities for persisting and accessing
genus-level genome dataset artifacts in the Ampire data hierarchy.
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import Dict, Tuple

# Third-Party Library Imports
import pandas as pd

# App Imports
from app.utils.config_utils import load_yaml, write_yaml
from app.utils.table_utils import load_table_with_columns


# Genomes Registry
def register_genus_genomes_dataset(
    *,
    genus_name: str,
    data_root: Path,
    assemblies: pd.DataFrame,
    parameters: Dict,
    overwrite: bool = False,
) -> Tuple[Path, Path]:
    """
    Register genome dataset artifacts for a given genus.

    This function persists genus-level genome download results
    into the Ampire data hierarchy, including assembly metadata
    and dataset provenance information.

    Parameters
    ----------
    genus_name : str
        Scientific genus name (e.g., "Escherichia").
    data_root : pathlib.Path
        Root directory for data storage (e.g., data/interim).
    assemblies : pd.DataFrame
        Assembly-level metadata table produced by genome
        retrieval pipelines (e.g., NCBI RefSeq downloads).
    parameters : dict
        Parameters used during genome retrieval, recorded
        for dataset reproducibility.
    overwrite : bool, optional
        Whether to overwrite an existing registered dataset.
        If False and the dataset already exists, an exception
        is raised.

    Returns
    -------
    Tuple[Path, Path]
        (dataset_yaml_path, assemblies_csv_path) registered
        under the Ampire data hierarchy.

    Raises
    ------
    FileExistsError
        If the genomes dataset has already been registered
        and overwrite is set to False.
    """

    # Construct paths
    genus_root = data_root / genus_name
    genomes_dir = genus_root / "genomes"
    dataset_yaml = genomes_dir / "dataset.yaml"
    assemblies_csv = genomes_dir / "assemblies.csv"

    # Check existing registration
    if dataset_yaml.exists() and assemblies_csv.exists() and not overwrite:
        raise FileExistsError(
            "Genomes dataset artifacts already exist and overwrite is disabled."
        )

    # Create directories (only now we touch filesystem)
    genomes_dir.mkdir(parents=True, exist_ok=True)

    # Write dataset metadata (YAML)
    dataset_manifest = {
        "dataset": {
            "name": "genomes",
            "scope": "genus",
            "genus": genus_name,
        },
        "source": {
            "tool": "ncbi-genome-download",
            "database": "NCBI RefSeq",
        },
        "parameters": parameters,
        "provenance": {
            "retrieved_at": pd.Timestamp.utcnow().isoformat(),
            "registered_by": "register_genus_genomes_dataset",
        },
        "artifacts": {
            "assemblies": "assemblies.csv",
            "fasta_dir": "fasta/",
        },
    }
    write_yaml(dataset_yaml, dataset_manifest, overwrite=overwrite)

    # Write assemblies table (CSV)
    assemblies.to_csv(
        assemblies_csv,
        index=False,
    )

    return dataset_yaml, assemblies_csv


def is_genus_genomes_registered(
    genus_name: str,
    data_root: Path,
) -> bool:
    """
    Determine whether genome dataset artifacts for a genus
    have been registered.

    Parameters
    ----------
    genus_name : str
        Scientific genus name.
    data_root : pathlib.Path
        Root directory for data storage (e.g., data/interim).

    Returns
    -------
    bool
        True if both genome dataset artifacts are present;
        False otherwise.
    """
    # Construct paths
    genomes_dir = data_root / genus_name / "genomes"
    assemblies_csv = genomes_dir / "assemblies.csv"
    dataset_yaml = genomes_dir / "dataset.yaml"

    return assemblies_csv.is_file() and dataset_yaml.is_file()


def load_genus_genomes_dataset(
    genus_name: str,
    data_root: Path,
) -> Tuple[Dict, pd.DataFrame]:
    """
    Load registered genome dataset artifacts for a genus.

    Parameters
    ----------
    genus_name : str
        Scientific genus name.
    data_root : pathlib.Path
        Root directory for data storage (e.g., data/interim).

    Returns
    -------
    Tuple[dict, pd.DataFrame]
        (dataset_metadata, assemblies) loaded from the
        registered genome dataset artifacts.

    Raises
    ------
    FileNotFoundError
        If the genome dataset has not been registered and
        required artifacts are missing.
    """
    # Construct paths
    genomes_dir = data_root / genus_name / "genomes"
    dataset_yaml = genomes_dir / "dataset.yaml"
    assemblies_csv = genomes_dir / "assemblies.csv"

    # Check existence
    if not dataset_yaml.is_file() or not assemblies_csv.is_file():
        raise FileNotFoundError(f"Genomes dataset not registered: '{genus_name}'")

    # Load artifacts
    dataset_metadata = load_yaml(dataset_yaml)
    assemblies = load_table_with_columns(str(assemblies_csv))

    return dataset_metadata, assemblies
