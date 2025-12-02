#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

# Check input arguments
if [ $# -lt 2 ]; then
  echo "Usage: $0 <PROJECT_DIR> <MICROBIOME_DIR>"
  exit 1
fi

# Set paths
PROJECT_DIR=$1
MICROBIOME_DIR=$2
METADATA_INPUT="$MICROBIOME_DIR/hmp2_metadata_2018-08-20.csv"
METADATA_OUTPUT="$PROJECT_DIR/data/processed/HMP2/nfcore_metadata.tsv"
MICROBIOME_INPUT="$MICROBIOME_DIR/biopsy_16S/samplesheet.csv"
MICROBIOME_SILVA_TRAIN="$PROJECT_DIR/data/external/silva/nr99_v138.1_wSpecies_train_set.fa.gz"
MICROBIOME_SILVA_SPECIES="$PROJECT_DIR/data/external/silva/species_assignment_v138.1.fa.gz"
MICROBIOME_OUTPUT_DIR="$PROJECT_DIR/experiments/HMP2/biopsy_16S"
REL_ABUNDANCE_TABLE="$PROJECT_DIR/experiments/HMP2/biopsy_16S/qiime2/rel_abundance_tables/rel-table-6.tsv"
GENUS_OUTPUT="$PROJECT_DIR/data/processed/HMP2/genus_list.csv"
BACDIVE_CONFIG="$PROJECT_DIR/configs/bacdive.json"
BACDIVE_OUTPUT="$PROJECT_DIR/data/processed/HMP2/genus_bacdive_gram.csv"
LOG_FILE="$PROJECT_DIR/logs/microbiome.log"

# nfcore parameter
DADA2_THREADS=8
FASTP_THREADS=8
CUTADAPT_THREADS=8
CPUS=24
MEM="256.GB"
TIME="48.h"

# Create logs folder if missing
mkdir -p "$PROJECT_DIR/logs"

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
CONDA_ENV="Ampire"
conda activate $CONDA_ENV

# ===================== Step 1: Prepare metadata =====================
rm -f "$LOG_FILE"
rm -f "$METADATA_OUTPUT"
python "$PROJECT_DIR/src/main.py" \
  --stage prepare_hmp2_16s_metadata \
  --log_path "$LOG_FILE" \
  --hmp2_metadata_input_csv "$METADATA_INPUT" \
  --hmp2_metadata_output_tsv "$METADATA_OUTPUT"

# ===================== Step 2: Microbiome composition analysis =====================
rm -rf "$PROJECT_DIR/.nextflow.log"
rm -rf "$PROJECT_DIR/experiments/HMP2"
rm -rf "$PROJECT_DIR/work"
rm -rf "$PROJECT_DIR/.nextflow"
python "$PROJECT_DIR/src/main.py" \
  --stage microbiome_composition \
  --log_path "$LOG_FILE" \
  --microbiome_input_csv "$MICROBIOME_INPUT" \
  --microbiome_metadata "$METADATA_OUTPUT" \
  --microbiome_silva_train "$MICROBIOME_SILVA_TRAIN" \
  --microbiome_silva_species "$MICROBIOME_SILVA_SPECIES" \
  --microbiome_output_dir "$MICROBIOME_OUTPUT_DIR" \
  --microbiome_dada2_threads "$DADA2_THREADS" \
  --microbiome_fastp_threads "$FASTP_THREADS" \
  --microbiome_cutadapt_threads "$CUTADAPT_THREADS" \
  --microbiome_max_cpus "$CPUS" \
  --microbiome_max_memory "$MEM" \
  --microbiome_max_time "$TIME"

# ===================== Step 3: Extract genus taxonomy =====================
python "$PROJECT_DIR/src/main.py" \
  --stage extract_hmp2_genus_taxonomy \
  --log_path "$LOG_FILE" \
  --hmp2_rel_abundance_tsv "$REL_ABUNDANCE_TABLE" \
  --hmp2_genus_output_csv "$GENUS_OUTPUT"

# ===================== Step 4: Query BacDive for Gram classification of each genus =====================
python "$PROJECT_DIR/src/main.py" \
  --stage query_bacdive_genus \
  --log_path "$LOG_FILE" \
  --bacdive_input_csv "$GENUS_OUTPUT" \
  --bacdive_config "$BACDIVE_CONFIG" \
  --bacdive_output_csv "$BACDIVE_OUTPUT"

# Deactivate the Conda environment
conda deactivate