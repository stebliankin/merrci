#!/bin/bash
#echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting Calculating PTR..."

# ----------------------------------- Import variables ----------------------------------------------------------------

#source ../config.cfg

# ---------------------------------------------------------------------------------------------------------------------

GENOME_LIST_STRING=""
if [ -f "$GENOMES_LIST" ]; then

    for NZ in `cat $GENOMES_LIST`; do
        GENOME_LIST_STRING=$GENOME_LIST_STRING" "$FASTA_INDIVIDUAL_DIR"/"$NZ".fna"
    done

    #echo "List of FASTA files: " $GENOME_LIST_STRING
    rm -f $PTR_TSV
    rm -f $PTR_PDF

    #echo bPTR -f$GENOME_LIST_STRING -s $SAM_FILE -o $PTR_TSV -plot $PTR_PDF -m "coverage" >> $LOG_FILE

    bPTR -f$GENOME_LIST_STRING -s $SAM_DIR"/"$SAM_FILE -o $PTR_TSV -plot $PTR_PDF -m "coverage"
    #bPTR -f$GENOME_LIST_STRING -s $SAM_DIR"/"$SAM_FILE -o $PTR_TSV -plot $PTR_PDF -m "gc_skew"

fi
#echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] End Calculating PTR..."