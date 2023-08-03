import pandas

data = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)

logi = (data['allele_type'] == 'unknown') & (data['allele_description'] != 'unknown') & (data['allele_description'] != '')

data[logi].to_csv('unknowns_with_description.tsv', sep='\t', index=False)