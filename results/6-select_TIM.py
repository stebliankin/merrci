
# The purpose of this script is to select TIM and AMP/MEM cohorts

import pandas as pd

def get_cohort(row):
    if row['r_Ticarcillin-Clavulanate']>0:
        return 'TIM'
    elif row['r_Ampicillin']>0 or row['r_Meropenem']:
        return 'AMP/MEM'
    else:
        return 'Other'

df = pd.read_csv('./analysis-out/4-Select_major_ABR/PTR_species_filtered_metadata_major_AMR_RPKM.csv')

df['Treatment'] = df.apply(lambda row: get_cohort(row), axis=1)

#print(len(df))

df = df[df['Treatment']!='Other']

print("Number of samples in TIM cohort: {}".format(len(df[df['Treatment']=='TIM'])))
print("Number of samples in AMP/MEM cohort: {}".format(len(df[df['Treatment']=='AMP/MEM'])))

ABR_genes = [x for x in df.columns if 'ARO' in x]

df = df[['Treatment'] + ABR_genes]

for gene in ABR_genes:
    if df[gene].mean() == 0:
        df = df.drop(gene, axis=1)

df.to_csv('analysis-out/4-Select_major_ABR/abr_tim.csv', index=False)