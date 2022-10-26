"""
format the allele list to have unique lines for each allele, with the PMIDs as a list in the last column:

> Input file
SPAC1006.08	unknown	etd1	etd1-1		unknown	PMID:15933715
SPAC1006.08	unknown	etd1	etd1-1		unknown	PMID:7845361

> Printed output
SPAC1006.08	unknown	etd1	etd1-1		unknown	PMID:15933715,PMID:7845361

"""

import sys


def main(input_file):
    with open(input_file, 'r') as ins:
        # Omit first line
        ins.readline()
        # Rename columns
        print('\t'.join(['systematic_id', 'allele_description', 'gene_name', 'allele_name', 'allele_synonym', 'allele_type', 'reference']))
        line = ins.readline()
        next_line = ins.readline()
        pmids = list()
        while (next_line):
            ls = line.strip().split('\t')
            pmids.append(ls[-1])
            next_ls = next_line.strip().split('\t')
            if ls[3] != next_ls[3]:
                line_out = ls[:]
                line_out[-1] = ','.join(pmids)
                print('\t'.join(line_out))
                pmids = list()

            line = next_line
            next_line = ins.readline()


if __name__ == "__main__":
    main(sys.argv[1])
