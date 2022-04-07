#!/bin/bash

CONFIG_PATH=$1
source $CONFIG_PATH

python3 ${SCRIPTS_DIR}/postProcessing/combine_partitions.py --RESULTS_DIR $RESULTS_SAMPLES \
                                                            --PTR_OUT_DIR $COMBINED_PARTITIONS_DIR

python3 ${SCRIPTS_DIR}/postProcessing/combine_samples.py --COMBINED_PARTITIONS_DIR $COMBINED_PARTITIONS_DIR \
                                                        --OUT_PTR $PROJECT_DIR"/merged_ptr.txt" \
                                                        --OUT_ABUNDANCE $PROJECT_DIR"/merged_abundance.txt"

