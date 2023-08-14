import pandas
import re

file_name = 'data/sgd/alleles_raw.tsv'

data = pandas.read_csv('data/sgd/alleles_sgd_raw.tsv', sep='\t', na_filter=False)

regex_gene_name = '[a-zA-Z][a-zA-Z][a-zA-Z]\d+'

# remove the deletion alleles
deletion_alleles = data['allele_name'].str.contains(f'^{regex_gene_name}-Δ$', regex=True)
data = data[~deletion_alleles].copy()

# Extract the descriptions that are before the ;
data['description_semicolon'] = data['allele_description'].str.split(';').str[0]


def extract_description_from_allele_name(allele_name):
    # Alleles that are ase1-123 should not be accounted
    if re.match(f'^{regex_gene_name}-\d+$', allele_name):
        return ''
    # ts alleles neither
    if re.match(f'^{regex_gene_name}-ts$', allele_name):
        return ''

    if re.match(f'^{regex_gene_name}-', allele_name):
        return re.sub(f'^{regex_gene_name}-', '', allele_name)
    return ''


data['description_name'] = data['allele_name'].apply(extract_description_from_allele_name)
exclude_allele_type = data['allele_type'].str.contains('frameshift')
data.loc[~exclude_allele_type, 'allele_type'] = 'amino_acid_dummy'

# Merge solutions from different PMIDs that are the same
groupby_columns = list(data.columns)
groupby_columns.remove('reference')
data = data.groupby(groupby_columns, as_index=False).agg({'reference': ','.join})


alleles_description_name = data[data['description_name'] != ''].copy()
alleles_description_name.drop(columns=['allele_description', 'description_semicolon'], inplace=True)
alleles_description_name.rename(columns={'description_name': 'allele_description'}, inplace=True)
alleles_description_name.to_csv('data/sgd/alleles_description_name.tsv', sep='\t', index=False)

alleles_description_semicolon = data[data['description_semicolon'] != ''].copy()
alleles_description_semicolon.drop(columns=['allele_description', 'description_name'], inplace=True)
alleles_description_semicolon.rename(columns={'description_semicolon': 'allele_description'}, inplace=True)
alleles_description_semicolon.to_csv('data/sgd/alleles_description_semicolon.tsv', sep='\t', index=False)
