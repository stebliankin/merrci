#!/bin/bash
# Objective:
#   Run split fasta
#
# Author: Vitalii Stebliankin (vsteb002@fiu.edu)

source ${SCRIPTS_DIR_LOCAL}/preprocessing/configLocal.sh


FASTA=$KRAKEN_DB_LOCAL"/library/bacteria/library-clean.fna"
OUT=$DB_PARTITIONS_LOCAL

mkdir -p $OUT

python ${SCRIPTS_DIR_LOCAL}/preprocessing/preprocessingModules/split_fasta_file.py --file $FASTA \
                            --partitions "64" \
                            --affix "kraken" \
                            --version "92" \
                            --out $OUT