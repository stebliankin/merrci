import pandas as pd
from scipy.stats import pearsonr


#### Constants
klebsiella_name_merrci = 'Klebsiella pneumoniae#abundance'
klebsiella_name_metaphlan = 's__Klebsiella_pneumoniae'

ecoli_name_merrci = 'Escherichia coli#abundance'
ecoli_name_metaphlan = 's__Escherichia_coli'

staph_name_merrci = 'Staphylococcus epidermidis#abundance'
staph_name_metaphlan = 's__Staphylococcus_epidermidis'

merrci_abundance = './data/PTR_species_filtered_metadata.csv'
metaphlan_abundance = './out_abundance/merged_abundance_subset.csv'

merrci_df = pd.read_csv(merrci_abundance, sep=',')
merrci_df = merrci_df[['sample', ecoli_name_merrci, klebsiella_name_merrci, staph_name_merrci]]
metaphlan_df = pd.read_csv(metaphlan_abundance, sep=',')
metaphlan_df[klebsiella_name_metaphlan] = metaphlan_df[klebsiella_name_metaphlan] /100
metaphlan_df[ecoli_name_metaphlan] = metaphlan_df[ecoli_name_metaphlan] /100
metaphlan_df[staph_name_metaphlan] = metaphlan_df[staph_name_metaphlan] /100

merged_df = merrci_df.merge(metaphlan_df, how='left', on='sample')


#### Compute Correlation
corr_klebsiella, pval = pearsonr(list(merged_df[klebsiella_name_merrci]), list(merged_df[klebsiella_name_metaphlan]))
print("Pearson correlation for Klebsiella pneumoniae: {}; p-value: {}".format(corr_klebsiella, pval))


#### Compute Correlation
corr_ecoli, pval = pearsonr(list(merged_df[ecoli_name_merrci]), list(merged_df[ecoli_name_metaphlan]))
print("Pearson correlation for Escherichia coli: {}; p-value: {}".format(corr_ecoli, pval))

#### Compute Correlation
corr_staph, pval = pearsonr(list(merged_df[staph_name_merrci]), list(merged_df[staph_name_metaphlan]))
print("Pearson correlation for Staphylococcus epidermidis: {}; p-value: {}".format(corr_staph, pval))


