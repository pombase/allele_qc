import pandas
import sys


def main(input_file):

    # Get the file name without extension
    input_file_name = input_file.split('.')[0]
    all_alleles = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)
    data = pandas.read_csv(input_file, sep='\t', na_filter=False)

    # merge on systematic_id and allele_name
    data = pandas.merge(data, all_alleles[['systematic_id', 'allele_name', 'allele_description', 'gene_name', 'allele_synonym']], on=['systematic_id', 'allele_name', 'allele_description'], how='left')

    for c1, c2 in [('change_description_to', 'allele_description'), ('change_name_to', 'allele_name'), ('change_type_to', 'allele_type')]:
        data.loc[data[c1] != '', c2] = data.loc[data[c1] != '', c1]
        data.drop(c1, axis=1, inplace=True)

    data.to_csv(input_file_name + '_formatted_for_check.tsv', sep='\t', index=False)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        main(arg)