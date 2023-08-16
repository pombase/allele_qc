"""
Filters out the alleles without phenotype annotations based on the canto allele list.
Just to make sure that there is no allele in canto that has no annotations, but has an annotation in a phaf file,
we also use the phaf files.
It overwrites the 1st argument alleles.tsv file.

usage: python alleles.tsv alleles_canto.tsv alleles_phaf.tsv
"""

import pandas
import sys


def main(alleles_file, canto_alleles_file, phaf_alleles_file):
    alleles = pandas.read_csv(alleles_file, sep='\t', na_filter=False)

    canto_alleles = pandas.read_csv(canto_alleles_file, sep='\t', na_filter=False)
    canto_alleles.rename(columns={'gene_systematic_id': 'systematic_id'}, inplace=True)

    phaf_alleles = pandas.read_csv(phaf_alleles_file, sep='\t', na_filter=False)
    phaf_alleles.columns = ['systematic_id', 'allele_description', 'gene_name', 'allele_name', 'allele_synonym', 'allele_type', 'reference']

    for df in (alleles, canto_alleles, phaf_alleles):
        df['unique_id'] = df['systematic_id'] + '$' + df['allele_name'] + '$' + df['allele_description']

    alleles_without_annotations = canto_alleles[(canto_alleles['annotation_count'] == 0) & ~canto_alleles['unique_id'].isin(set(phaf_alleles['unique_id']))]

    alleles = alleles[~alleles['unique_id'].isin(set(alleles_without_annotations['unique_id']))]
    alleles.to_csv(alleles_file, sep='\t', index=False)


if __name__ == '__main__':

    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
