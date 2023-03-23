import pandas

data = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)
exclude_names = list(pandas.read_csv('manual_alleles2remove.tsv', sep='\t', na_filter=False)['allele_name']) + \
    list(pandas.read_csv('manual_changes_formatted.tsv', sep='\t', na_filter=False)['allele_name']) + \
    list(pandas.read_csv('manual_fix_chimeras.tsv', sep='\t', na_filter=False)['allele_name'])

data = data.loc[(data.allele_type == 'other') & (~data.allele_name.isin(exclude_names)), ['systematic_id', 'allele_name', 'allele_description', 'reference']].copy()
data.sort_values(['reference', 'systematic_id', 'allele_name', 'allele_description'], inplace=True)

data['change_description_to'] = ''
data['change_name_to'] = ''
data['change_type_to'] = 'fusion_or_chimera'
data['comment'] = ''

data.to_csv('manual_fix_other.tsv', sep='\t', index=False)
