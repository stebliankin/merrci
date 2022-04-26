import pandas as pd

abr_genes = ['ARO.3000498ErmF',
'ARO.3003559cepA',
'ARO.3003318Streptomyces',
'ARO.3003730Bifidobacteria',
'ARO.3001214mdtM',
'ARO.3002868dfrG',
'ARO.3003318Streptomyces',
'ARO.3003730Bifidobacteria',
'ARO.3002062CMY-51',
'ARO.3002069CMY-59',
'ARO.3002079CMY-66',
'ARO.3002080CMY-67',
'ARO.3003773Clostridium',
'ARO.3000179tet(L)',
'ARO.3000195tetB(P)',
'ARO.3002175MIR-10',
'ARO.3004042Enterobacter',
'ARO.3002804FosA2',
'ARO.3003171ACT-36',
'ARO.3004042Enterobacter',
'ARO.3000616mel',
'ARO.3002438OKP-B-5',
'ARO.3002875dfrE',
'ARO.3003551emeA',
'ARO.3003949efrB',
"ARO:3004089ANT(3'')-IIa",
"ARO:3002556AAC(6')-Ii",
'ARO.3002819msrC',
'ARO.3000191tetQ',
'ARO.3000832evgA',
'ARO.3004039Escherichia',
'ARO.3002174MIR-9',
'ARO.3003548mdtN',
'ARO.3000024patA',
'ARO.3002985arnA',
'ARO.3003209FosA5',
'ARO.3000823ramA',
'ARO.3002399OXY-2-4',
'ARO.3001327mdtK',
'ARO.3002412OXY-5-2',
'ARO.3002415OXY-6-3',
"ARO.3002641APH(3')-Ia",
'ARO.3003056smeE',
'ARO.3004041Klebsiella',
'ARO.3004580Klebsiella',
'ARO.3004583Klebsiella',
'ARO.3003922oqxA',
'ARO.3001140SHV-86',
'ARO.3000190tetO',
'ARO.3000216acrB',
'ARO.3000815mgrA',
"ARO.3003905ANT(4')-Ib",
'ARO.3000617mecA',
'ARO.3002865dfrC',
'ARO.3000316mphA',
'ARO.3002405OXY-2-10',
'ARO.3003774Streptococcus']

df = pd.read_csv('analysis-out/4-Select_major_ABR/PTR_species_filtered_metadata_major_AMR_RPKM.csv')
abr_genes = ['ARO:'+x[4:] for x in abr_genes]
df = df[abr_genes]

df.to_csv('analysis-out/4-Select_major_ABR/selected_genes_causal.csv', index=False)