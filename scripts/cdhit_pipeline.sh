#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

# -------------------- Check input arguments --------------------
if [ $# -lt 1 ]; then
  echo "Usage: $0 <PROJECT_DIR>"
  exit 1
fi

# Set paths
PROJECT_DIR=$1
LOG_FILE="$PROJECT_DIR/logs/cdhit_dedup.log"
INPUT_FASTA="$PROJECT_DIR/data/processed/AMP/raw_merged/amp_raw.fasta"
METADATA_CSV="$PROJECT_DIR/data/processed/AMP/raw_merged/amp_raw_metadata.csv"
OUTPUT_DIR="$PROJECT_DIR/data/processed/AMP/dedup_100"

# CD-HIT parameter
IDENTITY=1.0
WORD_SIZE=5
MIN_LENGTH=5
THREADS=16
MEMORY_MB=0

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
conda activate Ampire

# ===================== Step 1: Run CD-HIT clustering =====================
rm -f "$LOG_FILE"
python "$PROJECT_DIR/src/main.py" \
  --stage cdhit_dedup \
  --log_path "$LOG_FILE" \
  --cdhit_input_fasta "$INPUT_FASTA" \
  --cdhit_metadata_csv "$METADATA_CSV" \
  --cdhit_output_dir "$OUTPUT_DIR" \
  --identity "$IDENTITY" \
  --word_size "$WORD_SIZE" \
  --min_length "$MIN_LENGTH" \
  --threads "$THREADS" \
  --memory_mb "$MEMORY_MB"

# Deactivate the Conda environment
conda deactivate
