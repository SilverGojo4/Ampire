# Ampire

This is the official codebase for **Ampire**.

(Antimicrobial Peptide Atlas – Genome-centric Exploration)

&nbsp;
(Preprint)[![Preprint](https://img.shields.io/badge/preprint-in%20preparation-lightgrey)](#) &nbsp;
(Documentation)[![Docs](https://img.shields.io/badge/docs-in%20progress-lightgrey)](#) &nbsp;
(PyPI)[![PyPI version](https://badge.fury.io/py/ampire.svg)](https://pypi.org/project/ampire) &nbsp;

**[2026.01.01]** Initial research pipeline under active development.

## Research Motivation

A growing body of experimental evidence indicates that many antimicrobial peptides
(AMPs) originate from bacteria, either as independently encoded peptides,
small open reading frames (smORFs), or proteolytically processed fragments of
larger proteins.

Despite this, most existing AMP resources focus on peptide-level properties
and prediction, with limited attention to _where_ these sequences may be encoded
within bacterial genomes.

**Ampire** adopts a genome-centric perspective by systematically mapping
experimentally validated AMPs to large collections of bacterial genomes,
with an initial focus on gut-associated bacteria.
Rather than asserting expression or biological activity, Ampire treats
sequence similarity as evidence of _putative genomic encoding potential_.

## Project Scope (v0.1)

The current version of Ampire focuses on building a reproducible research
pipeline for large-scale AMP–genome sequence analysis, including:

- Curation and deduplication of experimentally validated AMP sequences
- Downloading and organizing bacterial genome datasets
- Large-scale translated sequence alignment between AMPs and genomes
- Aggregation and summarization of alignment results at genome and species levels

Ampire does **not** claim that detected genomic matches are necessarily expressed,
processed, or biologically active. All results should be interpreted as
_putative sequence evidence_ requiring further experimental validation.

Future extensions may include interactive exploration tools and user-submitted
sequence analysis, but the current repository is intentionally scoped as a
research and analysis engine.

## Installation

Ampire is developed and tested with **Python >= 3.12**.

We recommend using Conda for environment management, especially when working
with external bioinformatics tools.

```bash
conda create -n ampire python=3.12
conda activate ampire
pip install -e .
```

## Project Structure

- `configs/` — Pipeline and stage-level configuration files
- `data/`
  - `raw/` — Unmodified source data (genomes, AMP databases)
  - `interim/` — Intermediate artifacts (translations, indices)
  - `processed/` — Final analysis-ready results
- `src/app/` — Core pipeline implementation
  - `processing/` — Genome download, translation, and sequence alignment
  - `analysis/` — Result aggregation and statistics
  - `visualization/` — Figures and summary plots
- `notebooks/` — Exploratory analysis and reporting
- `experiments/` — Experimental analyses and parameter exploration

## Contributing

We welcome contributions related to reproducibility, scalability,
and biological interpretation of AMP–genome analysis pipelines.

Please open an issue to discuss proposed changes before submitting a pull request.

## Acknowledgements

We sincerely thank the authors and maintainers of the following open-source projects
and resources that make this work possible:

- Biopython
- NCBI RefSeq
- BLAST / DIAMOND
- cd-hit
- ETE3
