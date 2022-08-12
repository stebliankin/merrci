#!/bin/bash

if [ -z "$1" ]; then
    echo "$(date) Error: please, specify the sample"
    echo "Usage: metaphlane3-single.sh \$SAMPLE"

else


    SAMPLE="$1"

    echo "[$(date)] Start Metaphlan3 for ${SAMPLE}"
    module load bowtie2-2.3.4.1-gcc-4.8.5-3io37rb

    SAMPLE_PATH="./samples/"


    OUT_ABUNDANCE_DIR_COHORT='./out_abundance/individual_samples/'
    mkdir -p $OUT_ABUNDANCE_DIR_COHORT

    OUT_ABUNDANCE=$OUT_ABUNDANCE_DIR_COHORT"/"$SAMPLE".txt"
    rm -f $OUT_ABUNDANCE

    #echo "${METAPHLAN_APP}/metaphlan2.py $SAMPLE_PATH --input_type fastq --no_map --ignore_archaea --ignore_eukaryotes --ignore_viruses --nproc $THREADS --tax_lev 's' > $OUT_ABUNDANCE"

    metaphlan ${SAMPLE_PATH}${SAMPLE}_1.fastq,${SAMPLE_PATH}${SAMPLE}_2.fastq --index latest --input_type fastq --no_map --ignore_archaea --ignore_eukaryotes --nproc $THREADS --tax_lev 's' -o $OUT_ABUNDANCE

    echo "[$(date)] End Metaphlan3 for ${SAMPLE}"
fi