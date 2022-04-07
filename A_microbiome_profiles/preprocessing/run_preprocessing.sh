#!/bin/bash

source configLocal.sh
# Prepare the database
export SCRIPTS_DIR_LOCAL

./preprocessingModules/1-build_refseq_db.sh
./preprocessingModules/2-clean-fasta.sh
./preprocessingModules/3-run_split_fasta.sh
python ${SCRIPTS_DIR_LOCAL}/preprocessing/preprocessingModules/4-extract_individual_fasta.py --INPUT_FASTA $KRAKEN_DB_LOCAL"/library/bacteria/library-clean.fna" \
															--OUT_DIR $FASTA_INDIVIDUAL_DIR_LOCAL
./preprocessingModules/5-index_partitions.sh
