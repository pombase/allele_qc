"""
Runs the main analysis pipeline (see readme). The output is stored in the file indicated
in --output, in addition, two extra output files are created, for the default output value:

Only the subset of alleles that needs fixing.
results/allele_results_errors.tsv

Only the subset of alleles that needs fixing, only the columns 'allele_description' and 'change_description_to'.
results/allele_results_errors_summarised.tsv
"""

from models import SyntaxRule, AllowedTypes
from refinement_functions import check_allele_description
from grammar import allowed_types_dict, composed_types_dict, aminoacid_grammar, nucleotide_grammar, disruption_grammar
import pickle
import pandas
import argparse
from common_autofix_functions import print_warnings
from genome_functions import handle_systematic_id_for_allele_qc


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


def check_fun(row, genome, syntax_rules_aminoacids, syntax_rules_nucleotides, syntax_rules_disruption, allowed_types):

    systematic_id = handle_systematic_id_for_allele_qc(row['systematic_id'], row['allele_name'], genome)
    if systematic_id is None:
        return empty_dict() | {'needs_fixing': True, 'invalid_error': 'systematic_id not in genome'}

    gene = genome[systematic_id]

    if 'amino_acid' in row['allele_type'] or 'nonsense_mutation' == row['allele_type']:
        if 'peptide' not in gene:
            return empty_dict() | {'needs_fixing': True, 'invalid_error': 'peptide sequence missing'}
        else:
            return check_allele_description(row.allele_description, syntax_rules_aminoacids, row.allele_type, allowed_types, gene)
    elif 'nucleotide' in row['allele_type']:
        return check_allele_description(row.allele_description, syntax_rules_nucleotides, row.allele_type, allowed_types, gene)
    elif 'disruption' == row['allele_type']:
        # TODO: handle this better and refactor
        if row['allele_description'] != '':
            return check_allele_description(row.allele_description, syntax_rules_disruption, row.allele_type, allowed_types, gene)
        # Special case where the description  is empty
        else:
            out_dict = check_allele_description(row.allele_name, syntax_rules_disruption, row.allele_type, allowed_types, gene)
            # The name matches the pattern
            if out_dict['change_description_to'] != '':
                return out_dict
            else:
                return empty_dict()
    else:
        return empty_dict()


if __name__ == '__main__':
    class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
    parser.add_argument('--genome', default='data/genome.pickle', help='input: genome dictionary built from contig files.')
    parser.add_argument('--alleles', default='data/alleles.tsv', help='input allele dataset')
    parser.add_argument('--output', default='results/allele_results.tsv', help='output file, also creates two extra files with the extension _errors.tsv and _errors_summarised.tsv')
    args = parser.parse_args()

    with open(args.genome, 'rb') as ins:
        genome = pickle.load(ins)

    allele_data = pandas.read_csv(args.alleles, delimiter='\t', na_filter=False)
    syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
    syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]
    syntax_rules_disruption = [SyntaxRule.parse_obj(r) for r in disruption_grammar]
    allowed_types = AllowedTypes(allowed_types=allowed_types_dict, composed_types=composed_types_dict)

    extra_cols = allele_data.apply(lambda row: check_fun(row, genome, syntax_rules_aminoacids, syntax_rules_nucleotides, syntax_rules_disruption, allowed_types), axis=1, result_type='expand')
    output_data = pandas.concat([allele_data, extra_cols], axis=1)
    column_order = ['systematic_id', 'gene_name', 'allele_id', 'allele_name', 'allele_description', 'allele_type', 'reference', 'allele_parts', 'needs_fixing', 'change_description_to', 'rules_applied', 'pattern_error', 'invalid_error', 'sequence_error', 'change_type_to']
    output_data = output_data[column_order]

    print_warnings(output_data[(output_data['needs_fixing'] == True) & (output_data['pattern_error'] == '') & (output_data['allele_type'].str.contains('nucleot') | output_data['allele_type'].str.contains('amino'))])
    output_data.to_csv(args.output, sep='\t', index=False)

    root_output_name = args.output.split('.')[0]

    output_data[output_data['needs_fixing'] == True].to_csv(f'{root_output_name}_errors.tsv', sep='\t', index=False)
    output_data[output_data['needs_fixing'] == True][['allele_description', 'change_description_to']].to_csv(f'{root_output_name}_errors_summarised.tsv', sep='\t', index=False)
