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
from genome_functions import process_systematic_id
import re

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


def handle_systematic_id_for_qc(row, genome: dict) -> str:
    """
    Returns the right systematic_id for the allele:
    - For no multi-transcript (row.systematic_id in genome), return row.systematic_id
    - For multi-transcript in which the allele name starts with the primary name + .1, .2, etc, (e.g. zas1.2-V123A) return that transcript (SPBC1198.04c.2).
    - For other multi-transcript genes, return the first transcript, ending in .1 (SPBC1198.04c.1).
    """

    # If it's in the genome dictionary, return it
    if row['systematic_id'] in genome:
        return row['systematic_id']

    # Get the first multiple transcript id, if it is a multi-transcript gene
    try:
        first_multi_transcript = process_systematic_id(row['systematic_id'], genome, 'first')
    except ValueError:
        return None

    # If we have reached here, it means that the systematic_id is from a multi-transcript gene
    # If the allele name contains the primary name .1, .2, etc, (e.g. zas1.2) then we pick that transcript (SPBC1198.04c.2).
    # Otherwise, we pick the first transcript

    transcript_id_regex = '^' + row['gene_name'] + '\.(\d+)'
    match = re.search(transcript_id_regex, row['allele_name'])
    if match:
        return row['systematic_id'] + '.' + match.groups()[0]
    return first_multi_transcript


def check_fun(row, genome, syntax_rules_aminoacids, syntax_rules_nucleotides, syntax_rules_disruption, allowed_types):

    systematic_id = handle_systematic_id_for_qc(row, genome)
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

    extra_cols = allele_data.apply(lambda row: check_fun(row, genome, syntax_rules_aminoacids, syntax_rules_nucleotides, syntax_rules_disruption, allowed_types), axis=1, result_type='expand')
    output_data = pandas.concat([allele_data, extra_cols], axis=1)
    column_order = ['systematic_id', 'allele_description', 'gene_name', 'allele_name', 'allele_synonym', 'allele_type', 'reference', 'allele_parts', 'needs_fixing', 'change_description_to', 'rules_applied', 'pattern_error', 'invalid_error', 'sequence_error', 'change_type_to']
    output_data = output_data[column_order]

    print_warnings(output_data[(output_data['needs_fixing'] == True) & (output_data['pattern_error'] == '') & (output_data['allele_type'].str.contains('nucleot') | output_data['allele_type'].str.contains('amino'))])
    output_data.to_csv(args.output, sep='\t', index=False)

    root_output_name = args.output.split('.')[0]

    output_data[output_data['needs_fixing'] == True].to_csv(f'{root_output_name}_errors.tsv', sep='\t', index=False)
    output_data[output_data['needs_fixing'] == True][['allele_description', 'change_description_to']].to_csv(f'{root_output_name}_errors_summarised.tsv', sep='\t', index=False)
