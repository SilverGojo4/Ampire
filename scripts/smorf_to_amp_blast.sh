#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

# -------------------- Check input arguments --------------------
if [ $# -lt 2 ]; then
  echo "Usage: $0 <PROJECT_DIR> <GENOMES_ROOT>"
  exit 1
fi

# Set paths
PROJECT_DIR=$1
GENOMES_ROOT=$2
MASTER_LOG="$PROJECT_DIR/logs/smorf_to_amp_master.log"
GENUS_LIST="$PROJECT_DIR/data/processed/HMP2/genus_list.csv"
BLAST_DB="$PROJECT_DIR/data/external/amp_db/amp_nr100"

# CD-HIT parameter
THREADS=8
EVALUE=5
WORD_SIZE=2

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
    GENUS_PATH="$GENOMES_ROOT/$GENUS_FOLDER"
    GENUS_LOG_FILE="$GENUS_PATH/${GENUS_LOG_NAME}.log"

    # Check folder existence
    if [ ! -d "$GENUS_PATH" ]; then
        echo "[SKIP] $GENUS_NAME — directory not found" | tee -a "$MASTER_LOG"
        continue
    fi

    # Run pipeline
    if python "$PROJECT_DIR/src/main.py" \
        --stage smorf_to_amp \
        --log_path "$GENUS_LOG_FILE" \
        --genomes_genus_name "$GENUS_NAME" \
        --genomes_output_dir "$GENUS_PATH" \
        --blast_reference_db "$BLAST_DB" \
        --blast_threads "$THREADS" \
        --blast_evalue "$EVALUE" \
        --blast_word_size "$WORD_SIZE"; then

        echo "[ OK ] $GENUS_NAME" | tee -a "$MASTER_LOG"
    else
        echo "[FAIL] $GENUS_NAME" | tee -a "$MASTER_LOG"
    fi
done

# Deactivate the Conda environment
conda deactivate
