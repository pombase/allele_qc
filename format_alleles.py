"""
format the allele list to have unique lines for each allele, with the PMIDs as a list in the last column:

> Input file
SPAC1006.08	unknown	etd1	etd1-1		unknown	PMID:15933715
SPAC1006.08	unknown	etd1	etd1-1		unknown	PMID:7845361

> Printed output
SPAC1006.08	unknown	etd1	etd1-1		unknown	PMID:15933715,PMID:7845361

To use: python format_alleles.py input_file.tsv output_file.tsv
"""

import sys
import pandas


def main(input_file, output_file):
    data = pandas.read_csv(input_file, sep='\t', na_filter=False)

    # Rename columns
    data.columns = ['systematic_id', 'allele_description', 'gene_name', 'allele_name', 'allele_synonym', 'allele_type', 'reference']

    # Join with comma-separated
    unique_identifiers = ['systematic_id', 'allele_description', 'allele_name']
    data2merge = data[unique_identifiers + ['reference']].groupby(unique_identifiers, as_index=False).agg({'reference': ','.join})
    data = data.drop(columns=['reference']).merge(data2merge, on=unique_identifiers).drop_duplicates()
    data.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
