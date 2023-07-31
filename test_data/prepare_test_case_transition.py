import pandas

data = pandas.read_csv('../results/allele_results.tsv', sep='\t', na_filter=False)

logi = (data.change_description_to == '') & (data.pattern_error == '') & (data.invalid_error == '') & data.allele_type.str.contains('amino|nucleotide', regex=True)

cols = ['systematic_id', 'allele_description', 'gene_name', 'allele_name', 'allele_synonym', 'allele_type']

data.loc[logi, cols].to_csv('../test_data/transition_test.tsv', sep='\t', index=False)
