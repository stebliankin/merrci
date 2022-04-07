# Objective:
# Combine abundance, PTR, metadata, and AMR

import pandas as pd
import numpy as np
import os
import math

# Step 1 - Combine abundance, PTR and metadata in format that each row is a sample
ptr_combined_file = "path/to/PTR_species_filtered_metadata_major.csv" #ptr averaged by sample + metadata

amr_abundance_dir = "path/to/AMR-abundance"

out_amr = "path/to/ptr_amr.csv"

ptr_combined_df = pd.read_csv(ptr_combined_file)

samples_list = list(ptr_combined_df["sample"].unique())
amr_dict={"sample":[]}
all_genes=[]

# Get list of all genes:
for sample in samples_list:
    # Read abundance file

    f = os.path.join(amr_abundance_dir,sample+"-abundance.txt")
    if os.path.exists(f):
        with open(f, 'r') as re:
            re.readline()
            for line in re.readlines():
                line = line.strip("\n")
                row = line.split("\t")
                gene, abundance = row[0], row[3]
                if gene not in all_genes:
                    all_genes.append(gene)
# assign dictionary:
for gene in all_genes:
    amr_dict[gene] = []

# Get all values
for sample in samples_list:
    # Read abundance file
    amr_dict["sample"].append(sample)
    f = os.path.join(amr_abundance_dir,sample+"-abundance.txt")
    if os.path.exists(f):
        genes_found = []
        with open(f, 'r') as re:
            re.readline()
            for line in re.readlines():
                line = line.strip("\n")
                row = line.split("\t")
                gene, abundance = row[0], float(row[3])
                genes_found.append(gene)
                amr_dict[gene].append(abundance)
        for gene in all_genes:
            if gene not in genes_found:
                amr_dict[gene].append(0)
    else:
        for gene in all_genes:
            amr_dict[gene].append(np.nan)

amr_df = pd.DataFrame(amr_dict)
amr_df.index = amr_df["sample"]
amr_df = amr_df.loc[:, (amr_df != 0).any(axis=0)]
# Scale amr genes:
amr_df_scaled = amr_df
amr_df_scaled = amr_df_scaled.drop("sample", axis=1)

for col in amr_df_scaled.columns:
    amr_df_scaled[col] = amr_df_scaled[col].apply(lambda x: math.log2(x+1))

amr_df_scaled["sample"] = amr_df_scaled.index
merged_df = ptr_combined_df.merge(amr_df_scaled, on="sample", how="left")
merged_df.to_csv(out_amr, index=False)
