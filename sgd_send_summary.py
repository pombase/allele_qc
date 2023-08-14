import pandas

data1 = pandas.read_csv('results/sgd/allele_description_name_qc_errors.tsv', sep='\t', na_filter=False)

data2 = pandas.read_csv('results/sgd/allele_description_semicolon_qc_errors.tsv', sep='\t', na_filter=False)

data2 = data2[~data2['allele_name'].isin(set(data1['allele_name']))].copy()

data = pandas.concat([data1, data2], ignore_index=True)
data = data[['systematic_id', 'gene_name', 'allele_name', 'allele_description', 'change_description_to', 'sequence_error']].copy()

# Merge with the original descriptions

original_data = pandas.read_csv('data/sgd/alleles_sgd_raw.tsv', sep='\t', na_filter=False)
original_data = original_data[['allele_name', 'allele_description']].drop_duplicates()
original_data.rename(columns={'allele_description': 'original_description'}, inplace=True)

data = data.merge(original_data, on='allele_name', how='left')
data[data['sequence_error'] != ''].to_csv('results/sgd/allele_sequence_errors.tsv', sep='\t', index=False)