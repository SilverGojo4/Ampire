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
LOG_FILE="$PROJECT_DIR/logs/microbiome.log"

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
conda activate Ampire

# ===================== Step 1: Prepare metadata =====================
rm -f "$LOG_FILE"
rm -f "$PROJECT_DIR/data/processed/HMP2/nfcore_metadata.tsv"
python "$PROJECT_DIR/src/main.py" \
  --stage prepare_hmp2_16s_metadata \
  --log_path "$LOG_FILE" \
  --hmp2_metadata_input_csv "$PROJECT_DIR/data/raw/HMP2/hmp2_metadata_2018-08-20.csv" \
  --hmp2_metadata_output_tsv "$PROJECT_DIR/data/processed/HMP2/nfcore_metadata.tsv"

# ===================== Step 2: Microbiome composition analysis =====================
rm -rf "$PROJECT_DIR/.nextflow.log"
rm -rf "$PROJECT_DIR/experiments/HMP2"
rm -rf "$PROJECT_DIR/work"
rm -rf "$PROJECT_DIR/.nextflow"
python "$PROJECT_DIR/src/main.py" \
  --stage microbiome_composition \
  --log_path "$LOG_FILE" \
  --microbiome_input_csv "$PROJECT_DIR/data/raw/HMP2/biopsy_16S/samplesheet.csv" \
  --microbiome_metadata "$PROJECT_DIR/data/processed/HMP2/nfcore_metadata.tsv" \
  --microbiome_silva_train "$PROJECT_DIR/data/external/silva/nr99_v138.1_wSpecies_train_set.fa.gz" \
  --microbiome_silva_species "$PROJECT_DIR/data/external/silva/species_assignment_v138.1.fa.gz" \
  --microbiome_output_dir "$PROJECT_DIR/experiments/HMP2/biopsy_16S" \
  --microbiome_dada2_threads 8 \
  --microbiome_fastp_threads 8 \
  --microbiome_cutadapt_threads 8 \
  --microbiome_max_cpus 24 \
  --microbiome_max_memory "256.GB" \
  --microbiome_max_time "48.h"

# ===================== Step 3: Extract genus taxonomy =====================
python "$PROJECT_DIR/src/main.py" \
  --stage extract_hmp2_genus_taxonomy \
  --log_path "$LOG_FILE" \
  --hmp2_rel_abundance_tsv "$PROJECT_DIR/experiments/HMP2/biopsy_16S/qiime2/rel_abundance_tables/rel-table-6.tsv" \
  --hmp2_genus_output_csv "$PROJECT_DIR/data/processed/HMP2/genus_list.csv"

# ===================== Step 4: Query BacDive for Gram classification of each genus =====================
python "$PROJECT_DIR/src/main.py" \
  --stage query_bacdive_genus \
  --log_path "$LOG_FILE" \
  --bacdive_input_csv "$PROJECT_DIR/data/processed/HMP2/genus_list.csv" \
  --bacdive_config "$PROJECT_DIR/configs/bacdive.json" \
  --bacdive_output_csv "$PROJECT_DIR/data/processed/HMP2/genus_bacdive_gram.csv"

# Deactivate the Conda environment
conda deactivate