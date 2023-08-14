import pandas

data1 = pandas.read_csv('results/sgd/allele_description_name_qc_errors.tsv', sep='\t', na_filter=False)
data1['description_comes_from'] = 'name'

data1['unique_id'] = data1['allele_name'] + '$' + data1['allele_description']
data2 = pandas.read_csv('results/sgd/allele_description_semicolon_qc_errors.tsv', sep='\t', na_filter=False)
data2['unique_id'] = data2['allele_name'] + '$' + data2['allele_description']

data2 = data2[~data2['unique_id'].isin(set(data1['unique_id']))].copy()
data2['description_comes_from'] = 'description before semicolon'

data = pandas.concat([data1, data2], ignore_index=True)
data = data[['systematic_id', 'gene_name', 'allele_name', 'allele_description', 'change_description_to', 'sequence_error', 'description_comes_from']].copy()

# Merge with the original descriptions

original_data = pandas.read_csv('data/sgd/alleles_sgd_raw.tsv', sep='\t', na_filter=False)
original_data = original_data[['allele_name', 'allele_description']].drop_duplicates()
original_data.rename(columns={'allele_description': 'original_description'}, inplace=True)

data = data.merge(original_data, on='allele_name', how='left')
data.sort_values(by=['systematic_id', 'allele_name'], inplace=True)
data[data['sequence_error'] != ''].to_csv('results/sgd/allele_sequence_errors.tsv', sep='\t', index=False)