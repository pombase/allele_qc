"""
Runs the main analysis pipeline (see readme). The output is stored in the file indicated
in --output, in addition, two extra output files are created, for the default output value:

Only the subset of alleles that needs fixing.
results/allele_results_errors.tsv

Only the subset of alleles that needs fixing, only the columns 'allele_description' and 'change_description_to'.
results/allele_results_errors_summarised.tsv
"""

from models import SyntaxRule
from refinement_functions import check_allele_description
from grammar import allowed_types, aminoacid_grammar, nucleotide_grammar, disruption_grammar
import pickle
import pandas
import argparse
from common_autofix_functions import print_warnings


def empty_dict():
    """
    Return the dictionary with error
    """
    return {
        'allele_parts': '',
        'needs_fixing': False,
        'change_description_to': '',
        'rules_applied': '',
        'pattern_error': '',
        'invalid_error': '',
        'sequence_error': '',
        'change_type_to': ''
    }


def main(allele_data, genome, syntax_rules_aminoacids, syntax_rules_nucleotides, syntax_rules_disruption, allowed_types):
    output_data_list = list()
    for i, row in allele_data.iterrows():
        if 'amino_acid' in row['allele_type'] or 'nonsense_mutation' == row['allele_type']:
            if row.systematic_id not in genome or 'peptide' not in genome[row.systematic_id]:
                output_data_list.append(dict(row) | empty_dict() | {'needs_fixing': True, 'invalid_error': f'Peptide sequence of {row.systematic_id} missing (perhaps multiple transcripts)'})
            else:
                output_data_list.append(dict(row) | check_allele_description(row.allele_description, syntax_rules_aminoacids, row.allele_type, allowed_types, genome[row.systematic_id]))
        elif 'nucleotide' in row['allele_type']:
            if row.systematic_id not in genome:
                output_data_list.append(dict(row) | empty_dict() | {'needs_fixing': True, 'invalid_error': f'Nucleotide sequence of {row.systematic_id} missing'})
            else:
                output_data_list.append(dict(row) | check_allele_description(row.allele_description, syntax_rules_nucleotides, row.allele_type, allowed_types, genome[row.systematic_id]))
        elif 'disruption' == row['allele_type']:
            # TODO: handle this better and refactor
            if row.systematic_id not in genome:
                output_data_list.append(dict(row) | empty_dict() | {'needs_fixing': True, 'invalid_error': f'systematic_id {row.systematic_id} missing (perhaps multiple transcripts)'})
            elif row['allele_description'] != '':
                output_data_list.append(dict(row) | check_allele_description(row.allele_description, syntax_rules_disruption, row.allele_type, allowed_types, genome[row.systematic_id]))
            # Special case where the description  is empty
            else:
                out_dict = check_allele_description(row.allele_name, syntax_rules_disruption, row.allele_type, allowed_types, genome[row.systematic_id])
                # The name matches the pattern
                if out_dict['change_description_to'] != '':
                    output_data_list.append(dict(row) | out_dict)
                else:
                    output_data_list.append(dict(row) | empty_dict())
        else:
            output_data_list.append(dict(row) | empty_dict())

    return pandas.DataFrame.from_records(output_data_list)


if __name__ == '__main__':
    class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
    parser.add_argument('--genome', default='data/genome.pickle', help='genome dictionary built from contig files.')
    parser.add_argument('--alleles', default='data/alleles.tsv')
    parser.add_argument('--output', default='results/allele_results.tsv')
    args = parser.parse_args()

    with open(args.genome, 'rb') as ins:
        genome = pickle.load(ins)

    allele_data = pandas.read_csv(args.alleles, delimiter='\t', na_filter=False)
    syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
    syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]
    syntax_rules_disruption = [SyntaxRule.parse_obj(r) for r in disruption_grammar]

    output_data = main(allele_data, genome, syntax_rules_aminoacids, syntax_rules_nucleotides, syntax_rules_disruption, allowed_types)
    print_warnings(output_data[(output_data['needs_fixing'] == True) & (output_data['pattern_error'] == '') & (output_data['allele_type'].str.contains('nucleot') | output_data['allele_type'].str.contains('amino'))])
    output_data.to_csv(args.output, sep='\t', index=False)

    root_output_name = args.output.split('.')[0]

    output_data[output_data['needs_fixing'] == True].to_csv(f'{root_output_name}_errors.tsv', sep='\t', index=False)
    output_data[output_data['needs_fixing'] == True][['allele_description', 'change_description_to']].to_csv(f'{root_output_name}_errors_summarised.tsv', sep='\t', index=False)
