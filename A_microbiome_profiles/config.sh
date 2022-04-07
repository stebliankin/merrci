#!/bin/bash
# ---------------------------------------------------------------------------------------------------------------------
# Config file for RAPToR pipeline
#
# Author:
#	Vitalii Stebliankin (vsteb002@fiu.edu)
#           Florida International University
#           Bioinformatics Research Group
# ---------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------
# Path Setup - PLEASE MODIFY THIS SECTION
# ---------------------------------------------------------------------------------------------------------------------
# Path setup
SCRIPTS_DIR="/path/to/RAPToR"
TMP_SCRIPTS_DIR="/path/to/temporary/files"
mkdir -p $TMP_SCRIPTS_DIR

LOGS_DIR="/path/to/logs"
mkdir -p $LOGS_DIR

PROJECT_DIR="path/to/result/folder"
mkdir -p $PROJECT_DIR

RAPToR_DB="path/to/RAPToR/database"
mkdir -p $RAPToR_DB

# Samples:
SAMPLES_DIR="/path/to/samples"
M1_AFFIX="_1.fastq"
M2_AFFIX="_2.fastq"
SAMPLES_FILE="path/to/list/of/samples"

BOWTIE2_APP="path/to/bowtie2/app"

THREADS="8"

# Account setup
SLURM_ACC="Slurm.account.name"
SLURM_QOS="QOS.credentials"

# ---------------------------------------------------------------------------------------------------------------------
# Automatic directory setup
# ---------------------------------------------------------------------------------------------------------------------
FASTA_INDIVIDUAL_DIR=$RAPToR_DB"/fasta-individual"
mkdir -p $FASTA_INDIVIDUAL_DIR
DB_PARTITIONS=$RAPToR_DB"/partitions"
mkdir -p $DB_PARTITIONS

DB_PREFIX="kraken_v92"

RESULTS_SAMPLES=$PROJECT_DIR
COMBINED_PARTITIONS_DIR=$PROJECT_DIR"/combined_partitions"
mkdir -p $COMBINED_PARTITIONS_DIR

# ABR Profiles
FASTA_FILE_ABR="path/to/downloaded/fasta/nucleotide_fasta_protein_homolog_model.fasta" # Make sure to download this file before using the pipeline
INDEX_ABR="/path/to/bwt2_index/" #will be computed in "1-index_AMR"
OUT_ABR="path/to/AMR-abundance"
SAM_DIR_ABR="path/to/SAM/alignment/file"
mkdir -p $SAM_DIR_ABR
# ---------------------------------------------------------------------------------------------------------------------
# Tool parameters
# ---------------------------------------------------------------------------------------------------------------------
BWT2_ACTION="TRUE"
ABUNDANCE_ACTION="TRUE"
GROWTH_ACTION="TRUE"

COV_CUTOFF="0.8"
GR_COV="3"
GRID_COV="1"
REASSIGN="FALSE"
METHOD="bPTR"
MODE="very-fast-local"