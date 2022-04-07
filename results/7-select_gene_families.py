# select variables of significance and merge it with gene families
import pandas as pd

causal_graph = './B-out-RPKM/bn_correlation.csv'

df = pd.read_csv(causal_graph)
df = df[df['from'].str.contains('ARO')]
df['from'] = df['from'].apply(lambda x: x[:11].replace('.',':'))

meta_df = pd.read_csv('metadata/card-data/aro_index.tsv', sep='\t')

df = df.merge(meta_df, how='left', left_on = 'from', right_on='ARO Accession')

meta_df = pd.read_csv('metadata/card-data/aro_categories_index.tsv', sep='\t')

df = df.merge(meta_df, how='left', on='Protein Accession')

df.to_csv('./B-out-RPKM/gene_families.csv', index=False)