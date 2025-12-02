#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

# Check input arguments
if [ $# -lt 2 ]; then
  echo "Usage: $0 <PROJECT_DIR> <PROTEOMES_ROOT>"
  exit 1
fi

# Set paths
PROJECT_DIR=$1
PROTEOMES_ROOT=$2
MASTER_LOG="$PROJECT_DIR/logs/download_proteomes_master.log"
GENUS_LIST="$PROJECT_DIR/data/processed/HMP2/genus_list.csv"

# keep nothing (default)
KEEP_RAW_FLAG=""
KEEP_SPLIT_FLAG=""

# If you want: keep both
# KEEP_RAW_FLAG="--keep_raw"
# KEEP_SPLIT_FLAG="--keep_split"

# Create logs folder if missing
mkdir -p "$PROJECT_DIR/logs"

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
CONDA_ENV="Ampire"
conda activate $CONDA_ENV

# Main Loop
rm -f "$MASTER_LOG"
tail -n +2 "$GENUS_LIST" | while IFS= read -r GENUS_NAME; do
    [[ -z "$GENUS_NAME" ]] && continue

    # Convert Genus name to folder name
    GENUS_FOLDER=$(echo "$GENUS_NAME" | sed 's/ /_/g' | sed 's/-/_/g')

    # Lowercase version for log file name
    GENUS_LOG_NAME=$(echo "$GENUS_FOLDER" | tr '[:upper:]' '[:lower:]')
    GENUS_PATH="$PROTEOMES_ROOT/$GENUS_FOLDER"
    GENUS_LOG_FILE="$GENUS_PATH/${GENUS_LOG_NAME}.log"

    # Run pipeline
    mkdir -p "$GENUS_PATH"
    if python "$PROJECT_DIR/src/main.py" \
         --stage download_proteomes_by_genus \
        --log_path "$GENUS_LOG_FILE" \
        --proteomes_genus_name "$GENUS_NAME" \
        --proteomes_output_dir "$GENUS_PATH"\
        $KEEP_RAW_FLAG \
        $KEEP_SPLIT_FLAG; then

        echo "[ OK ] $GENUS_NAME" | tee -a "$MASTER_LOG"
    else
        echo "[FAIL] $GENUS_NAME" | tee -a "$MASTER_LOG"
    fi

    # Random delay
    SLEEP_TIME=$((RANDOM % 8 + 3))   # 3–10 seconds
    sleep $SLEEP_TIME

done

# Deactivate the Conda environment
conda deactivate
