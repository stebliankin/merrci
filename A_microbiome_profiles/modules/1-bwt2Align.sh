#!/bin/bash
#----------------------------------------------------------------------------------------
# Objective:
#	1) Perform a Bowtie2 alignment of a single partition
#		using Kraken database
#		Result SAM files will be used to calculate growth rate

# Author:
#	Vitalii Stebliankin (vsteb002@fiu.edu)
#           Florida International University
#           Bioinformatics Research Group
#----------------------------------------------------------------------------------------

#echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
#echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting bowtie2 alignment..."


mkdir -p $SAM_DIR
mkdir -p $RESULTS_DIR

rm -f $SAM_DIR"/"$SAM_FILE
mkdir -p $SAM_DIR

echo "  [" `date '+%m/%d/%y %H:%M:%S'` "] Starting bowtie2 alignment parition "$PARTITION
echo $BOWTIE2_APP
echo $SAM_DIR"/"$SAM_FILE
echo $INDEX

$BOWTIE2_APP"/"bowtie2 --threads $THREADS --local --$MODE --no-discordant --no-mixed --no-unal -q --reorder \
	 -x $INDEX -1 $M1 -2 $M2 > $SAM_DIR"/"$SAM_FILE
echo "  [" `date '+%m/%d/%y %H:%M:%S'` "] Done with bowtie2 alignment parition "$PARTITION

# Single end mode:
#echo "  [" `date '+%m/%d/%y %H:%M:%S'` "] Starting bowtie2 alignment parition "$PARTITION
#bowtie2 --threads $THREADS --local --$MODE --no-discordant --no-mixed --no-unal -q --reorder \
#	 -x $INDEX -U $M1  > $SAM_DIR"/"$SAM_FILE
#echo "  [" `date '+%m/%d/%y %H:%M:%S'` "] Done with bowtie2 alignment parition "$PARTITION


#echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
#echo "[" `date '+%m/%d/%y %H:%M:%S'` "] End bowtie2 alignment..."
