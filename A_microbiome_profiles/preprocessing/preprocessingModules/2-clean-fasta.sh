#!/bin/bash
# Objective:
#   Clean FASTA header names
#
# Author: Vitalii Stebliankin (vsteb002@fiu.edu)

source ${SCRIPTS_DIR_LOCAL}/preprocessing/configLocal.sh

FASTA=$KRAKEN_DB_LOCAL"/library/bacteria/library.fna"
OUT=$KRAKEN_DB_LOCAL"/library/bacteria/library-clean.fna"

python ${SCRIPTS_DIR_LOCAL}/preprocessing/preprocessingModules/clean-fasta.py --file $FASTA \
                                                                                    --out $OUT