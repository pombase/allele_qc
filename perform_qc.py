from models import SyntaxRule
from refinement_functions import find_allele_parts
from grammar import allowed_types, grammar
import pickle
import sys
import pandas

# TODO needs fixing
# Build a dictionary PMID - curs
pmid2curs_dict = dict()
with open('data/pubs_and_session_ids.csv') as ins:
    for line in ins:
        pmid, curs = line.strip().split(',')
        pmid2curs_dict[pmid] = curs

with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)


def main(input_file: str):
    dict_list = list()
    syntax_rules = [SyntaxRule.parse_obj(r) for r in grammar]
    with open(input_file) as ins:
        ins.readline()
        for line in ins:
            systematic_id, allele_description, gene_name, allele_name, allele_synonym, allele_type, pmid = line.strip().split('\t')
            base_dict = {
                'systematic_id': systematic_id,
                'gene_name': gene_name,
                'allele_name': allele_name,
                'allele_type': allele_type,
                'allele_description': allele_description
            }
            if 'nucleotide' in allele_type:
                continue
            if systematic_id not in genome or 'translation' not in genome[systematic_id]:
                dict_list.append(
                    {
                        'allele_parts': '',
                        'needs_fixing': True,
                        'rename_to': '',
                        'rules_applied': '',
                        'pattern_error': '',
                        'invalid_error': 'several transcript or CDS missing',
                        'sequence_error': '',
                        'change_type_to': ''
                    } | base_dict
                )
            else:
                dict_list.append(base_dict | find_allele_parts(allele_description, syntax_rules, allele_type, allowed_types, genome[systematic_id]))

    pandas.DataFrame.from_records(dict_list).to_csv('results/allele_results.tsv', sep='\t')


if __name__ == "__main__":
    main(sys.argv[1])
