#!/bin/bash

conda create --name metaphlan2 python=3.7
source activate metaphlan2
conda install metaphlan2=2.7.8
conda install numpy
conda install biopython
conda install bowtie2