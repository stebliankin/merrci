#!/bin/bash
# Objective:
#	Download SRA samples

# Dependencies:
#	SRA toolkit
#	https://ncbi.github.io/sra-tools/fastq-dump.html

# Author: Vitalii Stebliankin (vsteb002@fiu.edu)

SAMPLES_DIR="./samples/"
mkdir -p $SAMPLES_DIR

SAMPLES_LIST="./data/SRA_list.txt"


for ID in $(cat $SAMPLES_LIST); do
    echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Downloading "$ID
	fastq-dump --split-files $ID --outdir $SAMPLES_DIR
done

#scp -r Antibiotics_Gibson2016 hpclogin01.fiu.edu:/home/vsteb002/MERRCI_additional/