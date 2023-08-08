import sys
import pandas
from refinement_functions import check_allele_description
from grammar import allowed_types_dict, aminoacid_grammar, nucleotide_grammar
from models import SyntaxRule
import pickle


def main(input_file):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
    syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]
    data = pandas.read_csv(input_file, sep='\t')
    data.fillna('', inplace=True)
    for i, line in data.iterrows():
        allele_type = line['change_type_to'] if line['change_type_to'] != '' else line['allele_type']
        allele_description = line['change_description_to'] if line['change_description_to'] != '' else line['allele_description']
        if 'amino' in allele_type:
            check = check_allele_description(allele_description, syntax_rules_aminoacids, allele_type, allowed_types_dict, genome[line['systematic_id']])
        else:
            check = check_allele_description(allele_description, syntax_rules_nucleotides, allele_type, allowed_types_dict, genome[line['systematic_id']])
        if check['needs_fixing']:
            print('\t'.join(line.tolist()))
            print('seq_error:', check['sequence_error'], 'change_description_to:', check['change_description_to'], 'change_type_to:', check['change_type_to'])
            print()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python check_manual_changes.py <file_name>")
        sys.exit(1)

    main(sys.argv[1])