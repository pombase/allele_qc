"""
format the allele list to have unique lines for each allele, with the PMIDs as a list in the last column:

> Input file
SPAC1006.08	unknown	etd1	etd1-1		unknown	PMID:15933715
SPAC1006.08	unknown	etd1	etd1-1		unknown	PMID:7845361

> Printed output
SPAC1006.08	unknown	etd1	etd1-1		unknown	PMID:15933715,PMID:7845361

To use: python format_alleles.py
"""

import pandas
import json


def main():

    with open('data/allele_type_mapping.json', 'r') as f:
        allele_type_mapping = json.load(f)

    data_canto = pandas.read_csv('data/alleles_pre_format_canto.tsv', sep='\t', na_filter=False)
    data_phaf = pandas.read_csv('data/alleles_pre_format_phaf.tsv', sep='\t', na_filter=False)

    # Format alleles from phaf
    data_phaf.columns = ['systematic_id', 'allele_description', 'gene_name', 'allele_name', 'allele_synonym', 'allele_type', 'reference']
    unique_identifiers = ['systematic_id', 'allele_description', 'allele_name']
    data2merge = data_phaf[unique_identifiers + ['reference']].groupby(unique_identifiers, as_index=False).agg({'reference': ','.join})
    data_phaf = data_phaf.drop(columns=['reference']).merge(data2merge, on=unique_identifiers).drop_duplicates()

    # Rename columns
    data_canto.rename(columns={'gene_systematic_id': 'systematic_id', 'references': 'reference', 'allele_synonyms': 'allele_synonym'}, inplace=True)

    # Map allele types
    data_canto['allele_type'] = data_canto['allele_type'].apply(lambda x: allele_type_mapping.get(x, x))
    data_canto = data_canto[(data_canto.allele_type != 'deletion') & (data_canto.allele_type != 'wild_type') & (data_canto.annotation_count > 0)].copy()

    # Include the rows from the canto file that do not exist in the phaf file (they have different values of systematic_id, description)
    # We prioritise the phaf file because it includes all synonyms, and addresses the case of missing names
    data_canto['identifier'] = data_canto['systematic_id'] + '$$$$' + data_canto['allele_description']
    data_phaf['identifier'] = data_phaf['systematic_id'] + '$$$$' + data_phaf['allele_description']
    data_canto = data_canto[~data_canto.identifier.isin(data_phaf.identifier)].copy()

    # Sort the columns, merge and save
    column_order = ['systematic_id', 'allele_description', 'gene_name', 'allele_name', 'allele_synonym', 'allele_type', 'reference']
    output_data = pandas.concat([data_canto[column_order], data_phaf[column_order]])
    output_data[column_order].sort_values(['systematic_id', 'allele_name', 'allele_description']).to_csv('data/alleles.tsv', sep='\t', index=False)


if __name__ == "__main__":
    main()
