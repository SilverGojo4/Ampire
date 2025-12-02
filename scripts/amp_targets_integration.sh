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
DBAMP_CSV="$PROJECT_DIR/data/processed/dbAMP/targets_normalized.csv"
DBAASP_CSV="$PROJECT_DIR/data/processed/DBAASP/targets_normalized.csv"
DRAMP_CSV="$PROJECT_DIR/data/processed/Dramp/targets_normalized.csv"
OUTPUT_DIR="$PROJECT_DIR/data/processed/AMP/raw_merged"
BACDIVE_CONFIG="$PROJECT_DIR/configs/bacdive.json"
BACDIVE_OUTPUT="$PROJECT_DIR/data/processed/AMP/raw_merged/targets_bacdive_gram.csv"
LOG_FILE="$PROJECT_DIR/logs/amp_targets.log"

# Create logs folder if missing
mkdir -p "$PROJECT_DIR/logs"

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
CONDA_ENV="Ampire"
conda activate $CONDA_ENV

# ===================== Step 1: Merge AMP targets across dbAMP / DBAASP / DRAMP =====================
rm -f "$LOG_FILE"
rm -f "$PROJECT_DIR/data/processed/merged/"
python "$PROJECT_DIR/src/main.py" \
  --stage merge_targets \
  --log_path "$LOG_FILE" \
  --merge_dbamp_csv "$DBAMP_CSV" \
  --merge_dbaasp_csv "$DBAASP_CSV" \
  --merge_dramp_csv "$DRAMP_CSV" \
  --merge_output_dir "$OUTPUT_DIR"

# ===================== Step 2: Annotate merged targets with BacDive metadata =====================
python "$PROJECT_DIR/src/main.py" \
  --stage query_bacdive_targets \
  --log_path "$LOG_FILE" \
  --bacdive_input_csv "$OUTPUT_DIR/targets.csv" \
  --bacdive_config "$BACDIVE_CONFIG" \
  --bacdive_output_csv "$BACDIVE_OUTPUT"

# Deactivate the Conda environment
conda deactivate