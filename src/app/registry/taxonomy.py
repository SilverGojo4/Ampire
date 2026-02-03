# pylint: disable=
"""
Registry utilities for persisting and accessing
genus-level taxonomy artifacts in the Ampire data hierarchy.
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


# Taxonomy Registry
def register_genus_taxonomy(
    *,
    genus_name: str,
    data_root: Path,
    genus_metadata: Dict,
    species_table: pd.DataFrame,
    overwrite: bool = False,
) -> Tuple[Path, Path]:
    """
    Register taxonomy artifacts for a given genus.

    This function persists pre-resolved taxonomy results into
    the Ampire data hierarchy.

    Parameters
    ----------
    genus_name : str
        Scientific genus name (e.g., "Escherichia").
    data_root : pathlib.Path
        Root directory for data storage (e.g., data/interim).
    genus_metadata : dict
        Genus-level taxonomy metadata returned by
        `resolve_species_from_genus`.
    species_table : pd.DataFrame
        Species-level taxonomy table returned by
        `resolve_species_from_genus`.
    overwrite : bool, optional
        Whether to overwrite existing registered taxonomy artifacts.
        If False and the taxonomy already exists, an exception is raised.

    Returns
    -------
    Tuple[Path, Path]
        (genus_yaml_path, species_csv_path) registered under
        the Ampire data hierarchy.

    Raises
    ------
    FileExistsError
        If taxonomy artifacts already exist and overwrite is False.
    """

    # Construct paths
    genus_root = data_root / genus_name
    taxonomy_dir = genus_root / "taxonomy"
    genus_yaml_path = taxonomy_dir / "genus.yaml"
    species_csv_path = taxonomy_dir / "species.csv"

    # Check existing registration
    if genus_yaml_path.exists() and species_csv_path.exists() and not overwrite:
        raise FileExistsError(
            "Genus taxonomy artifacts already exist and overwrite is disabled."
        )

    # Create directories (only now we touch filesystem)
    taxonomy_dir.mkdir(parents=True, exist_ok=True)

    # Write genus metadata (YAML)
    write_yaml(genus_yaml_path, genus_metadata, overwrite=overwrite)

    # Write species table (CSV)
    species_table.to_csv(
        species_csv_path,
        index=False,
    )

    return genus_yaml_path, species_csv_path


def is_genus_registered(
    genus_name: str,
    data_root: Path,
) -> bool:
    """
    Determine whether taxonomy artifacts for a genus are registered.

    Parameters
    ----------
    genus_name : str
        Scientific genus name.
    data_root : pathlib.Path
        Root directory for data storage (e.g., data/interim).

    Returns
    -------
    bool
        True if both genus-level and species-level taxonomy
        artifacts are present; False otherwise.
    """
    taxonomy_dir = data_root / genus_name / "taxonomy"
    genus_yaml = taxonomy_dir / "genus.yaml"
    species_csv = taxonomy_dir / "species.csv"

    return genus_yaml.is_file() and species_csv.is_file()


def load_genus_taxonomy(
    genus_name: str,
    data_root: Path,
) -> Tuple[Dict, pd.DataFrame]:
    """
    Load registered taxonomy artifacts for a genus.

    Parameters
    ----------
    genus_name : str
        Scientific genus name.
    data_root : pathlib.Path
        Root directory for data storage (e.g., data/interim).

    Returns
    -------
    Tuple[dict, pd.DataFrame]
        (genus_metadata, species_table) loaded from the
        registered taxonomy artifacts.

    Raises
    ------
    FileNotFoundError
        If the genus taxonomy has not been registered and
        required artifacts are missing.
    """
    # Construct paths
    taxonomy_dir = data_root / genus_name / "taxonomy"
    genus_yaml = taxonomy_dir / "genus.yaml"
    species_csv = taxonomy_dir / "species.csv"

    # Check existence
    if not genus_yaml.is_file() or not species_csv.is_file():
        raise FileNotFoundError(f"Genus taxonomy not registered: '{genus_name}'")

    # Load artifacts
    genus_metadata = load_yaml(genus_yaml)
    species_table = load_table_with_columns(
        str(species_csv),
        required_columns=["species_name", "species_taxid"],
    )

    return genus_metadata, species_table
