#!/bin/bash
#----------------------------------------------------------------------------------------
# Objective:
#	1) Download the RefSeq v92 database of compete reference genomes from Kraken


# Author:
#	Vitalii Stebliankin (vsteb002@fiu.edu)
#           Florida International University
#           Bioinformatics Research Group
#----------------------------------------------------------------------------------------

source ${SCRIPTS_DIR_LOCAL}/preprocessing/configLocal.sh

#kraken-build --download-library bacteria --db $KRAKEN_DB_LOCAL
kraken2-build --threads 24 --download-library bacteria --no-masking --db $KRAKEN_DB_LOCAL --skip-maps