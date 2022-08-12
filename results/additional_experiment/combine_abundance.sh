#!/bin/bash
# Objective:
#   Combine abundance files from all samples into a single abundance matrix

OUT_ABUNDANCE_DIR_COHORT=./out_abundance/individual_samples/

OUT=./out_abundance/merged_abundance.txt
rm -f $OUT

merge_metaphlan_tables.py ${OUT_ABUNDANCE_DIR_COHORT}/*.txt > $OUT
