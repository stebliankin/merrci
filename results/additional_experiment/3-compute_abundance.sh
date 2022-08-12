#!/bin/bash

source activate metaphlan3

SLURM_ACC=iacc_giri
SLURM_QOS=pq_giri

mkdir -p 'logs'
LOGS_METAPHLAN="./logs/metaphlan2"
mkdir -p $LOGS_METAPHLAN

ALL_METAPHLAN_JOBS=""
isFirstFlag=1
THREADS=4; export THREADS

while read SAMPLE; do
    METAPHLAN_SUBMIT=$(sbatch -J metaphlan \
            --account=$SLURM_ACC \
            --qos=$SLURM_QOS \
            -p investor \
            -N 1 \
            -n $THREADS \
            -o $LOGS_METAPHLAN"/stdout-${SAMPLE}.txt" \
            -e $LOGS_METAPHLAN"/stderr-${SAMPLE}.txt" \
            ./metaphlan3_single.sh $SAMPLE)
    echo $METAPHLAN_SUBMIT
    METAPHLAN_JID=$(echo $METAPHLAN_SUBMIT | awk -F' ' '{print $4}')
    if [ $isFirstFlag -eq 1 ]; then
        ALL_METAPHLAN_JOBS=$METAPHLAN_JID
        isFirstFlag=0
    else
        ALL_METAPHLAN_JOBS=$ALL_METAPHLAN_JOBS":"$METAPHLAN_JID
    fi
done < 'data/sra_list.txt'

echo $ALL_METAPHLAN_JOBS

LOGS_METAPHLAN_COMBINE="./logs/metaphlan2_combine"
mkdir -p $LOGS_METAPHLAN_COMBINE
# Combine abundance files
MET_COMBINE_SUBMIT=$(sbatch -J met-combine \
        --account=$SLURM_ACC \
        --qos=$SLURM_QOS \
        -p investor \
        -N 1 \
        -n $THREADS \
        -o $LOGS_METAPHLAN_COMBINE"/stdout-combine_abundance.txt" \
        -e $LOGS_METAPHLAN_COMBINE"/stderr-combine_abundance.txt" \
        --dependency=afterok:${ALL_METAPHLAN_JOBS} \
       ./combine_abundance.sh)
echo $MET_COMBINE_SUBMIT
MET_COMBINE_JID=$(echo $MET_COMBINE_SUBMIT | awk -F' ' '{print $4}')
echo $MET_COMBINE_JID