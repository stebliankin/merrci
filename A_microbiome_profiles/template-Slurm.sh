#!/bin/bash
#SBATCH -J SlurmPTR

# Number of nodes
#SBATCH -N 1

#SBATCH -p centos7

# Number of tasks
#SBATCH -n 16
#SBATCH -o LOGS_DIR_PLACEHOLDER/stdout-P-metaPTR_%a.txt
#SBATCH -e LOGS_DIR_PLACEHOLDER/stderr-P-metaPTR_%a.txt

# Array of jobs:
#SBATCH --array=1-PARTITIONS_PLACEHOLDER

#----------------------------------------------------------------------------------------
# Objective:
#	Submit jobs to calculate PTR for each partition
#
# Requirements:
#	* Python3
#	* bowtie2
#	* iRep
#
# Author:
#	Vitalii Stebliankin (vsteb002@fiu.edu)
#           Florida International University
#           Bioinformatics Research Group
#----------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------
#                                       Parameters setup
# ---------------------------------------------------------------------------------------------------------------------
echo "Running partition ${SLURM_ARRAY_TASK_ID}"

CONFIG="CONFIG_PLACEHOLDER"
source $CONFIG


SAMPLE="SAMPLE_PLACEHOLDER"
export SAMPLE

# Permanent
SAM_DIR=$PROJECT_DIR"/alignments"
echo  $SAM_DIR
mkdir -p $SAM_DIR
SAM_DIR=$SAM_DIR"/"$SAMPLE
mkdir -p $SAM_DIR
#SAM_DIR=$SAM_DIR"/"$MODE
#mkdir -p $SAM_DIR


# Environmental specific:
#source $CONDA_PATH"/activate" py3

PARTITION=${SLURM_ARRAY_TASK_ID}
export PARTITION

RESULTS_DIR=$RESULTS_SAMPLES"/"$SAMPLE
mkdir -p $RESULTS_DIR

#RESULTS_DIR=$RESULTS_DIR"/"$MODE


INDEX=$DB_PARTITIONS"/"$PARTITION"/"$DB_PREFIX
M1=$SAMPLES_DIR"/"$SAMPLE$M1_AFFIX
M2=$SAMPLES_DIR"/"$SAMPLE$M2_AFFIX
SAM_FILE=$PARTITION"_"$SAMPLE".sam"

ABUNDANCE_DIR=$RESULTS_DIR"/abundance"
mkdir -p $ABUNDANCE_DIR
ABUNDANCE_FILE=$ABUNDANCE_DIR"/"$PARTITION"_abundance_"$SAMPLE".txt"


GENOMES_LIST=$RESULTS_DIR"/"$PARTITION"_GenomesList_"$SAMPLE".txt"

PTR_TSV_DIR=$RESULTS_DIR"/bPTR_TSV"
mkdir -p $PTR_TSV_DIR

PTR_PDF_DIR=$RESULTS_DIR"/bPTR_PDF"
mkdir -p $PTR_PDF_DIR

PTR_TSV=$PTR_TSV_DIR"/"$PARTITION"_bPTR_"$SAMPLE".tsv"
PTR_PDF=$PTR_PDF_DIR"/"$PARTITION"_bPTR_"$SAMPLE".pdf"


echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting metaPTR pipeline. Sample: "$SAMPLE" Partition: "$PARTITION
## ---------------------------------------------------------------------------------------------------------------------
## Step 1 - Align
## ---------------------------------------------------------------------------------------------------------------------
if [ $BWT2_ACTION != "FALSE" ]; then
    echo "Starting bowtie2"
    . $SCRIPTS_DIR"/modules/1-bwt2Align.sh"
fi


# ---------------------------------------------------------------------------------------------------------------------
# Step 2 - Calculate Abundance
# ---------------------------------------------------------------------------------------------------------------------

if [ $ABUNDANCE_ACTION != "FALSE" ]; then
    . $SCRIPTS_DIR"/modules/2-abundance.sh"
fi

# ---------------------------------------------------------------------------------------------------------------------
# Step 3 - Calculate Growth Rates
# ---------------------------------------------------------------------------------------------------------------------

if [ $GROWTH_ACTION != "FALSE" ]; then

    . $SCRIPTS_DIR"/modules/3-PTR.sh"

fi

rm -f $GENOMES_LIST
rm -r $SAM_DIR"/"$SAM_FILE

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done with metaPTR pipeline. Sample: "$SAMPLE" Partition: "$PARTITION

