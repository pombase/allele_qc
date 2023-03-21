import pandas

data = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)
data = data[data.allele_type == 'other'].copy()
data = data.loc[data.allele_name.str.contains('fus|chim') | data.allele_description.str.contains('fus|chim'), ['systematic_id', 'allele_name', 'allele_description']]
data.to_csv('manual_fix_chimeras.tsv', sep='\t', index=False)
