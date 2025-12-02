#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

# Check input arguments
if [ $# -lt 1 ]; then
  echo "Usage: $0 <PROJECT_DIR>"
  exit 1
fi

# Set paths
PROJECT_DIR=$1
AMP_FASTA="$PROJECT_DIR/data/processed/AMP/dedup_100/amp_nr100.fasta"
DB_DIR="$PROJECT_DIR/data/external/amp_db"
OUTPUT_DB="$DB_DIR/amp_nr100"

# Ensure output directory exists
mkdir -p "$DB_DIR"

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
CONDA_ENV="Ampire"
conda activate $CONDA_ENV

# Build BLAST AMP DB
makeblastdb \
  -in "$AMP_FASTA" \
  -dbtype prot \
  -out "$OUTPUT_DB" \
  -parse_seqids

# Deactivate the Conda environment
conda deactivate
