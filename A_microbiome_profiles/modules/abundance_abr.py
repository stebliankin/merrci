#----------------------------------------------------------------------------------------
# Objective:
#	Calculate average coverage from alignment SAM file
#
# Input:
#   SAM file
#
# Output:
#   file with the following columns:
#       genome | coverage

# Author:
#	Vitalii Stebliankin (vsteb002@fiu.edu)
#           Florida International University
#           Bioinformatics Research Group
#----------------------------------------------------------------------------------------

import pandas as pd
import os
import argparse

#----------------------------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------------------------

def average_coverage_unit(genome_length, read_length):
    # Add average coverage based on single hit
    # C = NL / G
    #     where:
    #     C is the average coverage.
    #     N is the total number of reads that align to "this" sequence.
    #     L is the length of the read.
    #     G is the sequence length.
    return (int(read_length))/int(genome_length)



#----------------------------------------------------------------------------------------
# Initialization
#----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("--SAM",  help="Directory with SAM file")
parser.add_argument("--COV_CUTOFF", help="General coverage cutoff")
parser.add_argument("--N_READS", help="General coverage cutoff")
parser.add_argument("--OUTPUT_ABUNDANCE", help="Path to Output Abundance")

args = parser.parse_args()



SAM = args.SAM
COV_CUTOFF = float(args.COV_CUTOFF)
N_READS=int(args.N_READS)
OUTPUT_ABUNDANCE = args.OUTPUT_ABUNDANCE


if os.path.exists(OUTPUT_ABUNDANCE):
    os.remove(OUTPUT_ABUNDANCE)

# SAM = "../../data/control_1004.sam"
# READ_LENGTH = 250
# COV_CUTOFF = 0.8
# GR_COV = 5
# OUTPUT_ABUNDANCE="../../data/abundance.txt"
# OUTPUT_LIST="../../data/list.txt"
#


with open(SAM, 'r') as f:
    coverage_dict = {}
    length_dict = {}
    hits_per_million = {}
    for row in f.readlines():
        if ("@HD" in row) or ("@PG" in row):
            pass
        elif "@SQ" in row:
            row = row.strip("\n")
            row = row.split("\t")
            genome = row[1][3:]
            length = row[2].split(":")[1]
            length_dict[genome] = int(length)
            coverage_dict[genome] = 0
            hits_per_million[genome] = 0
        else:
            row = row.strip("\n")
            row = row.split("\t")
            genome = row[2]
            READ_LENGTH = len(row[9])
            curr_coverage = coverage_dict[genome] + average_coverage_unit(length_dict[genome], READ_LENGTH)
            coverage_dict[genome] = curr_coverage
            hits_per_million[genome] += 1000000 / N_READS

    coverage_df = pd.DataFrame.from_dict(coverage_dict, orient="index")
    coverage_df = coverage_df.rename(index=str, columns={0: "coverage"})
    coverage_df["genome"] = coverage_df.index
    coverage_df = coverage_df[coverage_df["coverage"]>COV_CUTOFF]
    coverage_df = coverage_df.reset_index(drop=True)

    hits_df = pd.DataFrame.from_dict(hits_per_million, orient="index")
    hits_df = hits_df.rename(index=str, columns={0: "hits_per_million"})
    hits_df["genome"] = hits_df.index
    hits_df = hits_df.reset_index(drop=True)

    coverage_df = coverage_df.merge(hits_df, how="left", on="genome")

    coverage_df = coverage_df.sort_values(by="coverage", ascending=False)
    summ = coverage_df["coverage"].sum()
    coverage_df["abundance"] = coverage_df["coverage"].apply(lambda x: x/summ)

    coverage_df[["genome", "coverage", "abundance", "hits_per_million"]].to_csv(OUTPUT_ABUNDANCE, sep="\t", index=False)
    # Getting the list of genomes that pass the coverage cutoff

