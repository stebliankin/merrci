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

rm -f $ABUNDANCE_FILE
rm -f $GENOMES_LIST

echo "  [" `date '+%m/%d/%y %H:%M:%S'` "] Starting calculating the abundance "$PARTITION



python3 $SCRIPTS_DIR"/utilities/abundance.py" --SAM $SAM_DIR"/"$SAM_FILE \
                                              --COV_CUTOFF $COV_CUTOFF --GR_COV $GR_COV \
                                              --OUTPUT_ABUNDANCE $ABUNDANCE_FILE \
                                              --OUTPUT_LIST $GENOMES_LIST \

echo "  [" `date '+%m/%d/%y %H:%M:%S'` "] Done with calculating abundance "$PARTITION




