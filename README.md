# Ampire

This is the official codebase for **Ampire**.

**Ampire** (Antimicrobial Peptide Atlas – Genome-centric Exploration) is a research
engine for large-scale, genome-centric analysis of antimicrobial peptides (AMPs).
It is designed for researchers and developers who aim to systematically explore
_putative genomic encoding potential_ of experimentally validated AMPs across
bacterial genomes.

&nbsp;
[![Preprint](https://img.shields.io/badge/preprint-in%20preparation-lightgrey)](#)
&nbsp;
[![Docs](https://img.shields.io/badge/docs-in%20progress-lightgrey)](#)
&nbsp;
[![PyPI version](https://badge.fury.io/py/ampire.svg)](https://pypi.org/project/ampire)

**Status:** _Initial research pipeline under active development_
**Last updated:** 2026-01-01

## Research Motivation

A growing body of experimental evidence indicates that many antimicrobial peptides
(AMPs) originate from bacteria, either as independently encoded peptides,
small open reading frames (smORFs), or proteolytically processed fragments of
larger proteins.

Despite this, most existing AMP resources and prediction frameworks focus on
peptide-level properties (e.g., physicochemical features or activity prediction),
with limited attention to **where these sequences may be encoded within bacterial
genomes**.

**Ampire** adopts a genome-centric perspective by systematically mapping
experimentally validated AMPs to large collections of bacterial genomes,
with an initial focus on gut-associated bacteria.
Rather than asserting expression or biological activity, Ampire treats
sequence similarity as evidence of **putative genomic encoding potential**.

This explicit distinction is intentional: Ampire is designed to support
hypothesis generation and large-scale exploratory analysis, not to replace
experimental validation.

## Project Scope (v0.1)

The current version of Ampire focuses on building a **reproducible research
pipeline** for large-scale AMP–genome sequence analysis, including:

- Curation and deduplication of experimentally validated AMP sequences
- Downloading, organizing, and registering bacterial genome datasets
- Large-scale translated sequence alignment between AMPs and genomes
- Aggregation and summarization of alignment results at genome and species levels

Ampire does **not** claim that detected genomic matches are necessarily expressed,
processed, or biologically active. All results should be interpreted as
_putative sequence evidence_ requiring further experimental validation.

Future extensions may include interactive exploration tools, dataset browsers,
and user-submitted sequence analysis. However, the current repository is
intentionally scoped as a **research and analysis engine**, serving as the
computational foundation for future platform development.

## Design Philosophy

Ampire is developed with the following principles in mind:

- **Genome-centric reasoning**
  Genomes, rather than peptides, are treated as the primary analytical context.

- **Reproducibility first**
  All analyses are driven by explicit configuration files, deterministic pipelines,
  and structured logging.

- **Clear epistemic boundaries**
  Sequence-level evidence is not conflated with expression, processing, or function.

- **Pipeline as reference, library as foundation**
  Core functionality is implemented as reusable modules, while pipelines serve as
  canonical reference implementations.

## Installation

Ampire is developed and tested with **Python >= 3.12**.

We recommend using Conda for environment management, especially when working
with external bioinformatics tools.

```bash
conda create -n ampire python=3.12
conda activate ampire
pip install -e .
```

## Configuration Overview

Ampire separates **public configuration interfaces** from **user-specific runtime
state**.

- `configs/examples/`
  Example configuration files that define the expected structure of inputs and
  pipeline stages. These files are version-controlled and safe to share.

- `configs/inputs/` _(user-specific)_
  Local input definitions (e.g., genus lists). These are intentionally excluded
  from version control.

- `configs/stages/` _(user-specific)_
  Pipeline execution configurations, including output locations and execution
  behavior.

- `configs/logging/`
  Centralized logging configuration defining Ampire’s logging behavior.

To get started:

```bash
cp configs/examples/inputs/genus.example.txt configs/inputs/genus.txt
cp configs/examples/stages/genomes_genus.example.yaml configs/stages/genomes_genus.yaml
```

Then edit the copied files according to your local environment.

## Current Pipelines

### `genomes.genus` — Genus-level genome dataset pipeline

The first canonical pipeline in Ampire builds reproducible, genus-level bacterial
genome datasets. It is responsible for:

- Resolving taxonomy and species membership for a given genus
- Downloading genome assemblies from public repositories
- Registering datasets and provenance metadata
- Producing standardized genome collections for downstream analysis

This pipeline serves as a **reference implementation** for future pipelines and
defines the standard structure for configuration, execution context, and logging
within Ampire.

## Project Structure

```
configs/        Pipeline configuration interfaces and examples
data/
  raw/          Unmodified source data (genomes, AMP databases)
  interim/      Intermediate artifacts (translations, indices)
  processed/    Final analysis-ready results
src/app/        Core pipeline and library implementation
  processing/   Genome download, translation, and alignment
  analysis/     Result aggregation and statistics
  visualization/Figures and summary plots
notebooks/      Exploratory analysis and reporting
experiments/    Experimental analyses and parameter exploration
```

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
