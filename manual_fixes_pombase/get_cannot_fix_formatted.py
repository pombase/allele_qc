import pandas

all_alleles = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)
data = pandas.read_csv('../results/allele_cannot_fix_sequence_errors.tsv', sep='\t', na_filter=False)

already_fixed = pandas.read_csv('manual_cannot_fix_new.tsv', sep='\t', na_filter=False)

data = data[~data['allele_name'].isin(already_fixed['allele_name'])].copy()

cannot_fix = pandas.read_csv('cannot_find.tsv', sep='\t', na_filter=False)
data = data[~data['allele_name'].isin(cannot_fix['allele_name'])].copy()

# merge on systematic_id and allele_name
data = pandas.merge(data, all_alleles, on=['systematic_id', 'allele_name', 'allele_description'], how='left')

data.sort_values(['reference', 'systematic_id', 'allele_name', 'allele_description'], inplace=True)

data['change_description_to'] = ''
data['change_name_to'] = ''
data['change_type_to'] = ''
data['comment'] = ''

columns = ['systematic_id', 'allele_name', 'allele_description', 'allele_type', 'reference', 'change_description_to', 'change_name_to', 'change_type_to', 'comment']

data[columns].to_csv('manual_cannot_fix.tsv', sep='\t', index=False)
