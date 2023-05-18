import pandas
import glob

data_list = list()

for f in glob.glob('allele_*.tsv'):
    d = pandas.read_csv(f, sep='\t', na_filter=False)
    if 'change_name_to' in d:
        data_list.append(d[['systematic_id', 'allele_name', 'change_name_to', 'change_type_to']])

data = pandas.concat(data_list)
data = data[(data.change_name_to != '') & (data.change_type_to != 'fusion_or_chimera')].copy()
data.rename(columns={'allele_name': 'Old name', 'change_name_to': 'New name'}, inplace=True)
data.drop(columns='change_type_to', inplace=True)
data.to_csv('temporary_changelog.tsv', sep='\t', index=False)

