"""
Extra analysis that takes the output file from 'perform_qc' (--alleles).

The analysis is applied to all alleles with references where sequence errors were found (--alleles_affected)
for genes where the coordinates have changed (--coordinate_changes). The output file contains two extra columns
when compared to --alleles:

    - 'after_coords_sequence_error': same meaning as 'sequence_error' column, but once the allele coordinates have been updated.
    - 'after_coords_change_description_to': same meaning as 'change_description_to' column, but once the allele coordinates have been updated.

"""

from models import SyntaxRule, find_rule
import argparse
import pandas
from genome_functions import get_other_index_from_alignment
from refinement_functions import build_regex2syntax_rule
import re
import json
from grammar import aminoacid_grammar
import pickle


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--alleles', default='results/allele_results.tsv')
parser.add_argument('--coordinate_changes', default='results/coordinate_changes_dict.json')
parser.add_argument('--alleles_affected', default='results/alleles_coordinate_change.tsv')
parser.add_argument('--output', default='results/allele_results_after_coordinates.tsv')
parser.add_argument('--genome', default='data/fasta_genome.pickle')
args = parser.parse_args()

with open(args.genome, 'rb') as ins:
    fasta_genome = pickle.load(ins)

with open(args.coordinate_changes, 'r') as ins:
    changes_dict = json.load(ins)

# Load alleles
allele_data = pandas.read_csv(args.alleles, delimiter='\t', na_filter=False)

# Exclude ambiguous
alleles_affected = pandas.read_csv(args.alleles_affected, delimiter='\t', na_filter=False)
allele_names_affected = set(alleles_affected['allele_name'][alleles_affected['uncertain_coordinate_change'] == False])

syntax_rules = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
regex2syntax_rule = build_regex2syntax_rule(syntax_rules)

allele_data['after_coords_sequence_error'] = ''
allele_data['after_coords_change_description_to'] = ''

for row_index, row in allele_data.iterrows():
    # Exclude
    if (row['allele_name'] not in allele_names_affected) or (row['systematic_id'] not in changes_dict):
        continue
    # If there was another error
    if not row['rules_applied']:
        continue
    rules_in_this_row = [find_rule(syntax_rules, rule_type, rule_name) for rule_type, rule_name in [rule_pair.split(':') for rule_pair in row['rules_applied'].split('|')]]
    new_alignment = changes_dict[row['systematic_id']]['new_alignment']
    old_alignment = changes_dict[row['systematic_id']]['old_alignment']

    new_allele_parts = list()
    new_sequence_errors = list()

    for rule, allele_part in zip(rules_in_this_row, row['allele_parts'].split('|')):
        if rule is None:
            print(row)
        match = re.match(rule.regex, allele_part)
        groups = list(match.groups())

        coordinate_error = list()
        for i in rule.coordinate_indexes:
            new_index = get_other_index_from_alignment(old_alignment, new_alignment, int(groups[i]))
            if new_index is None:
                # Special case where the N-term is truncated in a protein in which the first methionine was changed
                # making the protein shorter
                if groups[i] == '1' and 'partial' in row['allele_type']:
                    continue
                coordinate_error.append(groups[i])
                continue
            groups[i] = str(new_index)

        if len(coordinate_error):
            new_sequence_errors.append('coordinates do not exist: ' + ','.join(coordinate_error))
            continue

        new_sequence_errors.append(rule.check_sequence(groups, fasta_genome[row['systematic_id']]))
        new_allele_parts.append(rule.apply_syntax(groups))

    # If the change in coordinates does not change the allele name, skip
    new_allele_description = ','.join(new_allele_parts)
    if new_allele_description and ((new_allele_description == row['allele_description']) or (new_allele_description == row['change_description_to'])):
        continue

    allele_data.at[row_index, 'after_coords_sequence_error'] = '|'.join(new_sequence_errors) if any(new_sequence_errors) else ''
    allele_data.at[row_index, 'after_coords_change_description_to'] = new_allele_description
    allele_data.at[row_index, 'needs_fixing'] = True

allele_data.to_csv(args.output, sep='\t', index=False)

root_output_name = args.output.split('.')[0]

allele_data[allele_data['after_coords_change_description_to'] != ''].to_csv(f'{root_output_name}_concerned.tsv', sep='\t', index=False)
allele_data[allele_data['after_coords_change_description_to'] != ''][['allele_description', 'change_description_to', 'after_coords_change_description_to', 'sequence_error', 'after_coords_sequence_error']].to_csv(f'{root_output_name}_concerned_summarised.tsv', sep='\t', index=False)
