#-----------------------------------------------------------------------------------------------------------------------
# The purpose of this script is to combine PTR and abundance from all the samples
#
# Author:
#   Vitalii Stebliankin (vsteb002@fiu.edu)
#-----------------------------------------------------------------------------------------------------------------------

import pandas as pd
import argparse
import os

# PROJECT_DIR="/Users/stebliankin/Desktop/FLINTproject/data/JiInData"
# SAMPLES_LIST_FILE=PROJECT_DIR + "/WGS_Kidney_list.txt"

# PROJECT_DIR="/Users/stebliankin/Desktop/FLINTproject/data/antibiotics"
# SAMPLES_LIST_FILE=PROJECT_DIR + "/metadata/sra_list.txt"
#
# OUT_PTR = PROJECT_DIR+"/ptr_merged.txt"
# OUT_ABUNDANCE = PROJECT_DIR + "/abundance.txt"

parser = argparse.ArgumentParser()

parser.add_argument("--COMBINED_PARTITIONS_DIR", required=True, type=str, help="Path to merged PTR and abundance across paritions")
#parser.add_argument("--SAMPLES_LIST_FILE", required=True, type=str, help="List with samples")
parser.add_argument("--OUT_PTR", required=True, type=str, help="List with samples")
parser.add_argument("--OUT_ABUNDANCE", required=True, type=str, help="List with samples")

args = parser.parse_args()

COMBINED_PARTITIONS_DIR = args.COMBINED_PARTITIONS_DIR
OUT_PTR = args.OUT_PTR
OUT_ABUNDANCE = args.OUT_ABUNDANCE

def read_as_df(sample, COMBINED_PARTITIONS_DIR):
    sample = sample.strip("\n")
    ptr_file = COMBINED_PARTITIONS_DIR + "/bPTR_" + sample + ".txt"
    abundance_file = COMBINED_PARTITIONS_DIR + "/abundance_" + sample + ".txt"

    current_ptr_df = pd.read_csv(ptr_file, sep="\t", names=["genome", "ptr"])
    current_abundance_df = pd.read_csv(abundance_file, sep="\t", names=["genome", "coverage"])

    return current_ptr_df, current_abundance_df


def merge_time_series(COMBINED_PARTITIONS_DIR, OUT_PTR, OUT_ABUNDANCE):

    files = os.listdir(COMBINED_PARTITIONS_DIR)
    samples_list = []
    for f in files:
        sample = f.replace("abundance_","").replace("bPTR_","").replace(".txt","")
        if sample not in samples_list:
            samples_list.append(sample)
    for i, sample in enumerate(samples_list):
        sample = sample.replace("\n", "")
        #----------------------------
        # Read PTR and abundance
        #----------------------------

        current_ptr_df, current_abundance_df = read_as_df(sample, COMBINED_PARTITIONS_DIR)

        #------------------------------
        # Extract names for bPTR df
        #------------------------------
        # Extract NZ number:
            # From PTR
        current_ptr_df["NZ"] = current_ptr_df["genome"].apply(
            lambda row: row.split("/")[-1].split("_")[0] + "_" + row.split("/")[-1].split("_")[1])
        current_ptr_df = current_ptr_df[["NZ", "ptr"]]

            # From Abundance
        current_abundance_df["NZ"] = current_abundance_df["genome"].apply(lambda row: row.split("_")[1] + "_" + row.split("_")[2])

            # Get the full name from abundance file
        current_ptr_df = current_ptr_df.merge(current_abundance_df, how="left", on="NZ")

            # Rename columns
        current_ptr_df[sample] = current_ptr_df["ptr"]
        current_abundance_df[sample] = current_abundance_df["coverage"]

        # Normalize abundance
        sum_coverage = current_abundance_df[sample].sum()
        #current_abundance_df[sample] = current_abundance_df[sample].apply(lambda row: row/sum_coverage)



        current_ptr_df = current_ptr_df[["genome", sample]]
        current_abundance_df = current_abundance_df[["genome", sample]]

        if i==0:
            previous_ptr_df = current_ptr_df
            previous_abundance_df = current_abundance_df
        else:
            previous_ptr_df = previous_ptr_df.merge(current_ptr_df, on="genome", how="outer")
            previous_abundance_df = previous_abundance_df.merge(current_abundance_df, on="genome", how="outer")
        pass
        print(i)
        if i==18:
            pass

    previous_abundance_df.to_csv(OUT_ABUNDANCE, sep="\t", index=False)
    previous_ptr_df.to_csv(OUT_PTR, sep="\t", index=False)

#calculate_difference(PROJECT_DIR, SAMPLES_LIST_FILE, OUT_PTR, OUT_ABUNDANCE)

merge_time_series(COMBINED_PARTITIONS_DIR, OUT_PTR, OUT_ABUNDANCE)
        # Combine Abundance