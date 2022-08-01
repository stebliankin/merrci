#!/bin/bash
#----------------------------------------------------------------------------------------
# Objective:
#	Calculate abundance and coverage from alignment SAM file
#
# Input:
#   SAM file
#
# Output:
#   file with the following columns:
#       genome | abundance | coverage

# Author:
#	Vitalii Stebliankin (vsteb002@fiu.edu)
#           Florida International University
#           Bioinformatics Research Group
#----------------------------------------------------------------------------------------

#source ../config.cfg


mkdir -p $RESULTS_DIR

ID=$1

COV_CUTOFF=80

SAM_DIR="./SAM_files"
mkdir -p $SAM_DIR


echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Quantifying abundance for ${ID}..."

SAM_PATH="${SAM_DIR}/${ID}.sam"

## compute coverage
samtools sort $SAM_PATH > $SAM_PATH.sorted.sam
samtools coverage $SAM_PATH.sorted.sam > ${RESULTS_DIR}/${ID}-coverage.txt

let N_READS=$(wc -l $M1 | awk '{print $1}')/4

echo "  [" `date '+%m/%d/%y %H:%M:%S'` "] Number of reads: ${N_READS}"








