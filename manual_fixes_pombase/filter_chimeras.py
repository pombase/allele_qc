import pandas

data = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)
data = data[data.allele_type == 'other'].copy()
data = data.loc[data.allele_name.str.contains('fus|chim') | data.allele_description.str.contains('fus|chim'), ['systematic_id', 'allele_name', 'allele_description', 'reference']]
data.sort_values(['systematic_id', 'allele_name', 'allele_description', 'reference'], inplace=True)

data['change_description_to'] = ''
data['change_name_to'] = ''
data['change_type_to'] = 'fusion_or_chimera'
data['comment'] = ''

data.to_csv('manual_fix_chimeras.tsv', sep='\t', index=False)
