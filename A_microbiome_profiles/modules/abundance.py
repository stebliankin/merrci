#----------------------------------------------------------------------------------------
# Objective:
# Compute abundance
#
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

parser.add_argument("--ABUNDANCE_DIR",  help="Directory with pre-computed coverage")
parser.add_argument("--COV_CUTOFF", help="General coverage cutoff")
parser.add_argument("--N_READS", help="General coverage cutoff")
parser.add_argument("--OUTPUT_ABUNDANCE", help="Path to Output Abundance")

args = parser.parse_args()


abundance_dir = args.ABUNDANCE_DIR
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


combined_df = pd.DataFrame({"Samples#":[], "gene_name":[],
                            "gene_full_name":[],
                            "coverage":[], "numreads":[], "gene_length":[]})

for abundance_file in os.listdir(abundance_dir):
    if "-coverage.txt" in abundance_file:
        df_i = pd.read_csv(abundance_dir+abundance_file, sep='\t')
        df_i['gene_name'] = df_i['#rname'].apply(lambda x: x.split('|')[-1] if x.split('|')[-1] != 'RequiresSNPConfirmation' else x.split('|')[-2])

        df_i['gene_full_name'] = df_i['#rname']
        df_i = df_i.drop(['#rname'], axis=1)

        # remove 0 counts:
        df_i = df_i[df_i["numreads"]>0]

        # remove duplicates by selecting genes with maximum read counts:
        unique_df_i = df_i[df_i.groupby(['gene_name'])['coverage'].transform(max) == df_i['coverage']]

        # sometimes the read_count is exactly the same for several gene variants. In this case, randomly remove the duplicates:
        unique_df_i = unique_df_i.drop_duplicates(subset=['gene_name'])

        unique_df_i['Samples#'] = abundance_file.split('-coverage')[0]
        unique_df_i["gene_length"] = unique_df_i["endpos"] - unique_df_i["startpos"] + 1

        unique_df_i = unique_df_i[["Samples#", "gene_name",
                            "gene_full_name",
                            "coverage", "numreads", "gene_length"]]
        combined_df = combined_df.append(unique_df_i)

combined_df = combined_df.reset_index(drop=True)
#combined_df["Samples#"] = combined_df["Samples#"].apply(lambda x: x[:4])
combined_df["Function"] = combined_df['gene_full_name'].apply(lambda x: x.split('|')[2])

read_counts = pd.read_csv(N_READS)
read_counts["Samples#"] = read_counts["sample"]

combined_df = combined_df.merge(read_counts, how='left', on="Samples#")
# RPKM = ((Read_hits) * 10^3 * 10^6)/(Total_reads * gene_length)

combined_df["RPKM"] = (combined_df["numreads"] * 1000 * 1000000)/(combined_df["n_reads"]*combined_df['gene_length'])

def normalize_rpkm(combined_df):
    unique_samples = combined_df['Samples#'].unique()
    combined_df['sum'] = 0
    for sample in unique_samples:
        tmp_df = combined_df[combined_df['Samples#']==sample]
        combined_df['sum'] = combined_df.apply(lambda row: tmp_df['RPKM'].sum() if row['Samples#']==sample else row['sum'], axis=1)
    combined_df['abundance'] = combined_df.apply(lambda row: row['RPKM']*100/row['sum'], axis=1)
    return combined_df


combined_df = combined_df[combined_df["coverage"]>COV_CUTOFF]
combined_df = combined_df.fillna(0)
combined_df = normalize_rpkm(combined_df)

print("Total unique units identified: {}".format(len(combined_df['gene_name'].unique())))


combined_df.to_csv(OUTPUT_ABUNDANCE, index=False)

