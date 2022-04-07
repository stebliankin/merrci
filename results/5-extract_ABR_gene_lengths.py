from Bio import SeqIO
from collections import defaultdict
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')
#sns.set(rc={'figure.figsize':(11,2)})
#plt.figure(figsize=(5, 1))

ABR_FASTA_file = './metadata/card-data/nucleotide_fasta_protein_homolog_model.fasta'
LENGTH_OUT = './A-out/abr_gene_lengths.csv'

gene_lengths = defaultdict()

# Extract gene lenghts of ABR
with open(LENGTH_OUT, 'w') as out:
    out.write('ARO,Len\n')
    for record in SeqIO.parse(ABR_FASTA_file, "fasta"):
        aro_ID = record.description.split("|")[4]
        length = len(record.seq)
        gene_lengths[aro_ID] = length
        out.write('{},{}\n'.format(aro_ID, length))

# Extract list of ABR genes found in the dataset
#selected_ABR_file = 'analysis-out/4-Select_major_ABR/PTR_species_filtered_metadata_major_AMR.csv'
selected_ABR_file = 'A-out/old/ptr_amr.csv'

df = pd.read_csv(selected_ABR_file)
all_abr = [x.split('|')[4] for x in df.columns if 'ARO' in x]
all_length = [gene_lengths[x] for x in all_abr]

print(len(all_length))

# Plot histogram of gene lenghts
fig_dims = (8, 4)


# p = sns.displot(all_length, aspect=2, binwidth=50)
# p.fig.set_dpi(400)
# plt.savefig('analysis-out/figures/amr_genes_hist.png')

# Convert to RPKM
causal_input = './analysis-out/4-Select_major_ABR/PTR_species_filtered_metadata_major_AMR.csv'
converted_out = './analysis-out/4-Select_major_ABR/PTR_species_filtered_metadata_major_AMR_RPKM.csv'

df = pd.read_csv(causal_input)
abr_cols = [x for x in df.columns if 'ARO' in x]
for abr in abr_cols:
    # To convert CPM to RPKM:
    #   RPKM = (CPM * 10^3)/gene_length
    abr_len = gene_lengths[abr[:11]]
    df[abr] = df[abr].apply(lambda x: (x*1000)/abr_len)

df.to_csv(converted_out, index=False)


abr_file = 'A-out/old/ptr_amr.csv'
out_abr_file = 'A-out/ptr_amr_RPKM.csv'

df = pd.read_csv('A-out/old/ptr_amr.csv')
abr_cols = [x for x in df.columns if 'ARO' in x]
for abr in abr_cols:
    aro = abr.split('|')[4]
    abr_len = gene_lengths[aro]
    df[abr] = df[abr].apply(lambda x: (x*1000)/abr_len)

df.to_csv(out_abr_file, index=False)