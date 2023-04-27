import pandas

data = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)
data = data[data.allele_type == 'nonsense_mutation'].copy()

data.sort_values(['reference', 'systematic_id', 'allele_name', 'allele_description'], inplace=True)

data['change_description_to'] = ''
data['change_name_to'] = ''
data['change_type_to'] = 'partial_amino_acid_deletion'
data['comment_on_change'] = ''
data['add_synonym'] = ''
data['add_comment'] = ''


columns = ['systematic_id', 'allele_name', 'allele_description', 'allele_type', 'reference', 'change_description_to', 'change_name_to', 'change_type_to', 'comment_on_change', 'add_synonym', 'add_comment']

data[columns].to_csv('nonsense_mutations.tsv', sep='\t', index=False)
