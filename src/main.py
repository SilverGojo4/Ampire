# pylint: disable=line-too-long, import-error, wrong-import-position, broad-exception-caught, too-many-statements
"""
Ampire - Project Main Entry Point

This script serves as the centralized execution interface for the Ampire pipeline.
"""
# ============================== Standard Library Imports ==============================
import argparse
import importlib
import logging
import os
import sys

# ============================== Project Root Path Setup ==============================
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
LOGGING_PATH = os.path.join(BASE_PATH, "src/utils/logging_toolkit/src/python")

# Fix Python Path
if BASE_PATH not in sys.path:
    sys.path.append(BASE_PATH)
if LOGGING_PATH not in sys.path:
    sys.path.append(LOGGING_PATH)

# ============================== Project-Specific Imports ==============================
# Logging configuration and custom logger
from setup_logging import setup_logging

# ============================== Stage Configuration ==============================
SUPPORTED_STAGES = {
    "download_genomes_by_genus": {
        "title": "Download bacterial genomes from NCBI RefSeq by genus",
        "import_path": "src.preprocess.genomes.download_genomes.run_download_genomes_by_genus",
    },
    "smorf_by_genus": {
        "title": "Run smORFinder prediction pipeline by genus",
        "import_path": "src.preprocess.genomes.smorf_by_genus.run_smorf_by_genus",
    },
    "clean_dbamp_targets": {
        "title": "Clean dbAMP Targets into standardized long-format table",
        "import_path": "src.preprocess.amp.clean_dbamp_targets.run_clean_targets",
    },
    "clean_dbaasp_targets": {
        "title": "Clean DBAASP Targets into standardized long-format table",
        "import_path": "src.preprocess.amp.clean_dbaasp_targets.run_clean_targets",
    },
    "clean_dramp_targets": {
        "title": "Clean DRAMP Targets into standardized long-format table",
        "import_path": "src.preprocess.amp.clean_dramp_targets.run_clean_targets",
    },
    "normalize_targets": {
        "title": "Apply manually curated YAML mapping to dbAMP Targets",
        "import_path": "src.preprocess.amp.prepare_targets_mapping.run_normalize_targets",
    },
    "merge_targets": {
        "title": "Merge AMP targets across dbAMP, DBAASP, and DRAMP by sequence",
        "import_path": "src.preprocess.amp.merge_targets.run_merge_targets_by_sequence",
    },
    "prepare_hmp2_16s_metadata": {
        "title": "Extract HMP2 16S biopsy metadata for nf-core/ampliseq",
        "import_path": "src.preprocess.hmp2.prepare_metadata.run_extract_hmp2_16s_metadata",
    },
    "microbiome_composition": {
        "title": "Run nf-core/ampliseq for 16S microbiome composition analysis",
        "import_path": "src.analysis.hmp2.microbiome_composition.run_microbiome_composition",
    },
    "query_bacdive_targets": {
        "title": "Query BacDive for Gram classification of AMP targets",
        "import_path": "src.preprocess.amp.bacdive_gram_lookup.run_query_bacdive_targets",
    },
    "extract_hmp2_genus_taxonomy": {
        "title": "Extract genus-level taxonomy from HMP2 OTU table",
        "import_path": "src.preprocess.hmp2.extract_taxonomy.run_extract_hmp2_genus_taxonomy",
    },
    "query_bacdive_genus": {
        "title": "Query BacDive for genus-level Gram classification",
        "import_path": "src.preprocess.hmp2.bacdive_genus_lookup.run_query_bacdive_genus",
    },
    "cdhit_dedup": {
        "title": "Run CD-HIT for AMP deduplication",
        "import_path": "src.preprocess.amp.cdhit_dedup.run_cdhit_dedup",
    },
    "smorf_to_amp": {
        "title": "Run smORF-to-AMP BLAST screening by genus",
        "import_path": "src.analysis.amp_similarity.smorf_to_amp_blast.run_smorf_to_amp_by_genus",
    },
}


# ============================== Pipeline Dispatcher ==============================
def dispatch_stage(args: argparse.Namespace) -> None:
    """
    Dispatch execution to the appropriate pipeline stage using lazy import.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing stage-specific options.
    """
    # -------------------- Stage Validation --------------------
    # Ensure that the provided stage name is supported.
    stage = args.stage.lower()
    if stage not in SUPPORTED_STAGES:
        available = ", ".join(SUPPORTED_STAGES.keys())
        raise ValueError(f"Unknown stage '{stage}'. Available stages: {available}.")

    # -------------------- Log Path Validation --------------------
    # Ensure that the --log_path argument is provided.
    if not args.log_path:
        raise ValueError("Please specify '--log_path' to enable logging output.")

    # -------------------- Dynamic Stage Import --------------------
    # Perform a lazy import of the selected pipeline stage based on its import path.
    # This avoids loading all modules at startup and improves modularity.
    stage_info = SUPPORTED_STAGES[stage]
    module_path, func_name = stage_info["import_path"].rsplit(".", 1)
    module = importlib.import_module(module_path)
    stage_func = getattr(module, func_name)

    # -------------------- Logger Initialization --------------------
    # Set up the logging environment, loading the logging config and preparing output directory.
    log_config_file = os.path.join(BASE_PATH, "configs/logging.json")
    stage_log_path = os.path.abspath(args.log_path)
    stage_log_dir = os.path.dirname(stage_log_path)
    os.makedirs(name=stage_log_dir, exist_ok=True)

    logger = setup_logging(
        input_config_file=log_config_file,
        logger_name="combo_logger",
        handler_name="file",
        output_log_path=stage_log_path,
    )

    # -------------------- Log Pipeline Metadata --------------------
    logger.info("[ 'Pipeline Initialization Summary' ]")
    logger.log_pipeline_initialization(
        project_name="Ampire",
        line_width=120,
    )
    logger.add_spacer(level=logging.INFO, lines=1)

    # -------------------- Execute Stage Function --------------------
    # Dynamically call the selected stage function, passing all stage-specific arguments.
    extra_args = vars(args)
    extra_args.pop("stage", None)
    extra_args.pop("log_path", None)
    stage_func(base_path=BASE_PATH, logger=logger, **extra_args)
    logger.add_spacer(level=logging.INFO, lines=1)


# ============================== Main Entry ==============================
def main():
    """
    Main CLI entry point for the Ampire pipeline.
    Parses CLI arguments and routes execution to the selected pipeline stages.
    """
    # -------------------- Argument Parser --------------------
    available_stage_lines = [
        f"  - {stage:<15} {info['title']}" for stage, info in SUPPORTED_STAGES.items()
    ]
    available_stages_text = "\n".join(available_stage_lines)
    example_stage = list(SUPPORTED_STAGES.keys())[0]
    example_command = (
        f"  python main.py --stage {example_stage} --log_path logs/{example_stage}.log"
    )

    parser = argparse.ArgumentParser(
        description=(
            "Ampire - Accelerating the Discovery of Antimicrobial Peptides through Computational Intelligence Pipeline\n\n"
            "Available stages:\n"
            f"{available_stages_text}\n\n"
            "Example:\n"
            f"{example_command}\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # -------------------- General Options --------------------
    parser.add_argument(
        "--stage",
        type=str,
        choices=list(SUPPORTED_STAGES.keys()),
        required=True,
        help="Pipeline stage to run. Choose one of: "
        + ", ".join(SUPPORTED_STAGES.keys()),
    )
    parser.add_argument(
        "--log_path",
        type=str,
        required=True,
        help="Path for log file output.",
    )

    # -------------------- Options: genome download --------------------
    parser.add_argument(
        "--genomes_genus_name",
        type=str,
        help="Genus name to query (default: Escherichia).",
    )
    parser.add_argument(
        "--genomes_output_dir",
        type=str,
        help="Output dir for genomes (default: data/genomes/{genus_name}).",
    )

    # -------------------- Options: dbAMP target cleaning --------------------
    parser.add_argument(
        "--targets_input_path",
        type=str,
        help="Input dbAMP Excel/CSV (default: data/raw/dbAMP/dbAMP3_pepinfo.xlsx).",
    )
    parser.add_argument(
        "--targets_output_csv",
        type=str,
        help="Output cleaned CSV (default: data/processed/dbAMP/targets_clean.csv).",
    )

    # -------------------- Options: manual YAML mapping --------------------
    parser.add_argument(
        "--mapping_input_yaml",
        type=str,
        help="Input YAML mapping file (default: data/manual/targets_mapping/targets_mapping.yml).",
    )
    parser.add_argument(
        "--mapping_amp_input_csv",
        type=str,
        help="Input cleaned AMP targets CSV (default: data/processed/dbAMP/targets_clean.csv).",
    )

    # -------------------- Options: merge targets --------------------
    parser.add_argument(
        "--merge_dbamp_csv",
        type=str,
        help="Input dbAMP targets CSV (default: data/processed/dbAMP/targets_normalized.csv).",
    )
    parser.add_argument(
        "--merge_dbaasp_csv",
        type=str,
        help="Input DBAASP targets CSV (default: data/processed/DBAASP/targets_normalized.csv).",
    )
    parser.add_argument(
        "--merge_dramp_csv",
        type=str,
        help="Input DRAMP targets CSV (default: data/processed/DRAMP/targets_normalized.csv).",
    )
    parser.add_argument(
        "--merge_output_dir",
        type=str,
        help="Output directory for merged results (default: data/processed/merged/).",
    )

    # -------------------- Options: HMP2 16S metadata extraction --------------------
    parser.add_argument(
        "--hmp2_metadata_input_csv",
        type=str,
        help="Input HMP2 metadata CSV (default: data/raw/HMP2/hmp2_metadata_2018-08-20.csv).",
    )

    parser.add_argument(
        "--hmp2_metadata_output_tsv",
        type=str,
        help="Output TSV for nf-core/ampliseq (default: data/processed/HMP2/nfcore_metadata.tsv).",
    )

    # -------------------- Options: Microbiome composition (nf-core/ampliseq) --------------------
    parser.add_argument(
        "--microbiome_input_csv",
        type=str,
        help="Input sample sheet CSV for nf-core/ampliseq.",
    )

    parser.add_argument(
        "--microbiome_metadata",
        type=str,
        help="Sample metadata TSV for nf-core/ampliseq.",
    )

    parser.add_argument(
        "--microbiome_silva_train",
        type=str,
        help="Custom SILVA training set FASTA.",
    )

    parser.add_argument(
        "--microbiome_silva_species",
        type=str,
        help="Custom SILVA species FASTA.",
    )

    parser.add_argument(
        "--microbiome_output_dir",
        type=str,
        help="Output directory for nf-core/ampliseq results.",
    )

    parser.add_argument(
        "--microbiome_dada2_threads",
        type=int,
        help="Threads for DADA2 denoising.",
    )

    parser.add_argument(
        "--microbiome_fastp_threads",
        type=int,
        help="Threads for fastp preprocessing.",
    )

    parser.add_argument(
        "--microbiome_cutadapt_threads",
        type=int,
        help="Threads for cutadapt trimming.",
    )

    parser.add_argument(
        "--microbiome_max_cpus",
        type=int,
        help="Maximum CPU allocation.",
    )

    parser.add_argument(
        "--microbiome_max_memory",
        type=str,
        help="Max memory allocation (e.g., '64.GB').",
    )

    parser.add_argument(
        "--microbiome_max_time",
        type=str,
        help="Max time allocation (e.g., '48.h').",
    )

    # -------------------- Options: BacDive Gram classification --------------------
    parser.add_argument(
        "--bacdive_input_csv",
        type=str,
        help="Input CSV containing unique target names (default: data/processed/merged/targets_unique.csv).",
    )
    parser.add_argument(
        "--bacdive_config",
        type=str,
        help="Path to BacDive credential config JSON (default: configs/bacdive.json).",
    )
    parser.add_argument(
        "--bacdive_output_csv",
        type=str,
        help="Output CSV to save Gram classification results (default: data/processed/merged/targets_bacdive_gram.csv).",
    )

    # -------------------- Options: HMP2 genus taxonomy extraction --------------------
    parser.add_argument(
        "--hmp2_rel_abundance_tsv",
        type=str,
        help="Input HMP2 relative abundance TSV (default: data/processed/HMP2/rel_abundance_table.tsv).",
    )

    parser.add_argument(
        "--hmp2_genus_output_csv",
        type=str,
        help="Output CSV for genus list (default: data/processed/HMP2/genus_list.csv).",
    )

    # -------------------- Options: CD-HIT --------------------
    parser.add_argument(
        "--cdhit_input_fasta",
        type=str,
        help="Input FASTA for CD-HIT (default: data/processed/AMP/raw_merged/amp_raw.fasta).",
    )
    parser.add_argument(
        "--cdhit_metadata_csv",
        type=str,
        help="Input AMP metadata CSV (default: data/processed/AMP/raw_merged/amp_raw_metadata.csv).",
    )
    parser.add_argument(
        "--cdhit_output_dir",
        type=str,
        help="Output directory for CD-HIT results (default: data/processed/AMP/cdhit_nr100).",
    )
    parser.add_argument(
        "--identity",
        type=float,
        help="Identity threshold for CD-HIT (default: 1.0).",
    )
    parser.add_argument(
        "--word_size",
        type=int,
        help="CD-HIT word size (-n) (default: 5).",
    )
    parser.add_argument(
        "--min_length",
        type=int,
        help="Minimum sequence length threshold for CD-HIT (default: 5).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of threads for CD-HIT.",
    )
    parser.add_argument(
        "--memory_mb",
        type=int,
        help="Memory in MB for CD-HIT.",
    )

    # -------------------- Options: smORF-to-AMP BLAST screening --------------------
    parser.add_argument(
        "--blast_reference_db",
        type=str,
        help="Path prefix for AMP BLAST+ reference database (default: data/processed/AMP/blast_db/amp_db).",
    )

    parser.add_argument(
        "--blast_threads",
        type=int,
        help="Threads for BLASTP short-peptide search (default: 8).",
    )

    parser.add_argument(
        "--blast_evalue",
        type=float,
        help="E-value cutoff for BLASTP (default: 5).",
    )

    parser.add_argument(
        "--blast_word_size",
        type=int,
        help="Word size for short-peptide BLASTP (default: 2).",
    )

    args = parser.parse_args()

    # -------------------- Stage Execution --------------------
    try:
        dispatch_stage(args)
        print("Pipeline execution completed successfully.")
        sys.exit(0)

    except Exception as e:
        print(f"[Pipeline Error] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
