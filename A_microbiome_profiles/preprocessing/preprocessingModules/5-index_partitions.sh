#!/bin/bash
# ---------------------------------------------------------------------------------------------------------------------
#   OBJECTIVE:
#	The purpose of this script is to index a collection of bacterial genomes stored in a fasta file with the
#	bowtie2-build index builder program.
#
#   DEPENDENCIES:
#
#       • Bowtie2
#
#	AUTHORS:	Camilo Valdes (cvalde03@fiu.edu), Vitalii Stebliankin (vsteb002@fiu.edu)
#				Bioinformatics Research Group,
#				School of Computing and Information Sciences,
#				Florida International University (FIU)
#
# ---------------------------------------------------------------------------------------------------------------------
echo ""
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"
echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Starting..."
echo "[" `date '+%m/%d/%y %H:%M:%S'` "]"

source ${SCRIPTS_DIR_LOCAL}/preprocessing/configLocal.sh


NUMBER_OF_PARTITIONS=64
NUMBER_OF_THREADS=$THREADS

PARTITIONS_DIR=$DB_PARTITIONS_LOCAL

INDEX_NAME="kraken_v92"

PARTITIONS_ARRAY=( $(seq 1 $NUMBER_OF_PARTITIONS) )

for PARTITION in ${PARTITIONS_ARRAY[@]}; do
    echo "	[" `date '+%m/%d/%y %H:%M:%S'` "] indexing partition $PARTITION"
    OUT_DIR=$PARTITIONS_DIR"/$PARTITION"
    FASTA=$OUT_DIR"/"$INDEX_NAME".fasta"
    bowtie2-build -f --threads $NUMBER_OF_THREADS $FASTA $OUT_DIR"/"$INDEX_NAME
    rm $FASTA
done

echo "[" `date '+%m/%d/%y %H:%M:%S'` "] Done."
