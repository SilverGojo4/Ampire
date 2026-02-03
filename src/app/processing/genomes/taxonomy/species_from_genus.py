# pylint: disable=
"""
Taxonomy utilities for resolving species information from a genus name.

This module provides a stable, domain-level API for querying NCBI taxonomy
(via ETE3) and returning structured, reproducible taxonomy results.
"""

from __future__ import annotations

# Standard Library Imports
from datetime import datetime, timezone
from typing import Dict, List, Tuple

# Third-Party Library Imports
import pandas as pd
from ete3 import NCBITaxa


# Public API
def resolve_species_from_genus(
    genus_name: str,
) -> Tuple[Dict, pd.DataFrame]:
    """
    Resolve species-level taxonomy information from a genus name.

    Parameters
    ----------
    genus_name : str
        Scientific genus name (e.g., "Escherichia").

    Returns
    -------
    (genus_metadata, species_table) : Tuple[Dict, pd.DataFrame]
        genus_metadata : dict
            Genus-level taxonomy metadata with the following schema:
            {
              "genus": {
                "name": <str>,
                "taxid": <int>,
                "lineage": {
                  "ranks": [kingdom, phylum, class, order, family],
                  "names": [<str>, <str>, <str>, <str>, <str>]
                },
                "taxonomy_source": "NCBI Taxonomy (ETE3)",
                "retrieved_at": <ISO-8601 UTC timestamp>
              }
            }
        species_table : pd.DataFrame
            Species-level table with schema: species_name | species_taxid
            The table may be empty if the genus does not exist
            or has no species-level descendants.
    """
    ncbi = NCBITaxa()

    # Normalize input (taxonomy backend is case-sensitive in practice)
    genus_name = genus_name.strip()

    # Retrieve genus taxid
    name_map = ncbi.get_name_translator([genus_name])
    if not name_map or genus_name not in name_map:
        # Genus not found → return empty species table with minimal metadata
        genus_metadata = _build_empty_genus_metadata(genus_name)
        species_table = _empty_species_table()
        return genus_metadata, species_table

    genus_taxid = name_map[genus_name][0]

    # Build genus-level metadata
    genus_metadata = _build_genus_metadata(
        ncbi=ncbi,
        genus_name=genus_name,
        genus_taxid=genus_taxid,
    )

    # Resolve species under this genus
    species_table = _resolve_species_table(
        ncbi=ncbi,
        genus_taxid=genus_taxid,
    )

    return genus_metadata, species_table


# Internal Helpers
def _build_genus_metadata(
    ncbi: NCBITaxa,
    genus_name: str,
    genus_taxid: int,
) -> Dict:
    """
    Construct genus-level taxonomy metadata.

    Parameters
    ----------
    ncbi : NCBITaxa
        Initialized ETE3 taxonomy client.
    genus_name : str
        Scientific genus name.
    genus_taxid : int
        NCBI taxonomy identifier for the genus.

    Returns
    -------
    dict
        Genus-level taxonomy metadata dictionary.
    """
    lineage_taxids = ncbi.get_lineage(genus_taxid)
    rank_map = ncbi.get_rank(lineage_taxids)
    name_map = ncbi.get_taxid_translator(lineage_taxids)

    # Target lineage ranks (explicit contract)
    target_ranks: List[str] = ["kingdom", "phylum", "class", "order", "family"]

    lineage_ranks: List[str] = []
    lineage_names: List[str] = []

    for rank in target_ranks:
        for taxid in lineage_taxids:  # type: ignore
            if rank_map.get(taxid) == rank:
                lineage_ranks.append(rank)
                lineage_names.append(name_map.get(taxid))  # type: ignore
                break

    return {
        "genus": {
            "name": genus_name,
            "taxid": genus_taxid,
            "lineage": {
                "ranks": lineage_ranks,
                "names": lineage_names,
            },
            "taxonomy_source": "NCBI Taxonomy (ETE3)",
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
        }
    }


def _resolve_species_table(
    ncbi: NCBITaxa,
    genus_taxid: int,
) -> pd.DataFrame:
    """
    Resolve species-level descendants for a given genus taxid.

    Parameters
    ----------
    ncbi : NCBITaxa
        Initialized ETE3 taxonomy client.
    genus_taxid : int
        NCBI taxonomy identifier for the genus.

    Returns
    -------
    pd.DataFrame
        Species-level taxonomy table.
    """
    descendants = ncbi.get_descendant_taxa(
        genus_taxid,
        intermediate_nodes=True,
    )

    rank_map = ncbi.get_rank(descendants)
    name_map = ncbi.get_taxid_translator(descendants)

    records = [
        {
            "species_name": name_map[taxid],
            "species_taxid": taxid,
        }
        for taxid in descendants
        if rank_map.get(taxid) == "species"
    ]

    if not records:
        return _empty_species_table()

    return (
        pd.DataFrame.from_records(records)
        .sort_values("species_name")
        .reset_index(drop=True)
        .copy()
    )


def _empty_species_table() -> pd.DataFrame:
    """
    Return an empty species table with a stable schema.

    Returns
    -------
    pd.DataFrame
        Empty species-level DataFrame.
    """
    return pd.DataFrame(columns=["species_name", "species_taxid"])


def _build_empty_genus_metadata(
    genus_name: str,
) -> Dict:
    """
    Construct minimal genus metadata when the genus cannot be resolved.

    Parameters
    ----------
    genus_name : str
        Input genus name.

    Returns
    -------
    dict
        Minimal genus-level taxonomy metadata dictionary.
    """
    return {
        "genus": {
            "name": genus_name,
            "taxid": None,
            "lineage": {
                "ranks": [],
                "names": [],
            },
            "taxonomy_source": "NCBI Taxonomy (ETE3)",
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
        }
    }
