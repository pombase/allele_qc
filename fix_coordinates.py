from models import SyntaxRule, find_rule
import argparse
import pandas
from genome_functions import get_other_index_from_alignment
from refinement_functions import build_regex2syntax_rule
import re
import json
from grammar import aminoacid_grammar
from load_sequences import fasta_genome

parser = argparse.ArgumentParser(description='Build a dictionary of alignments based on the updated coordinates of genes')
parser.add_argument('--alleles', default='results/allele_results.tsv')
parser.add_argument('--coordinate_changes', default='results/coordinate_changes_dict.json')
parser.add_argument('--alleles_affected', default='results/alleles_coordinate_change.tsv')
parser.add_argument('--output', default='results/allele_results_after_coordinates.tsv')
args = parser.parse_args()


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
allele_data['after_coords_rename_to'] = ''

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

    allele_data.at[row_index, 'after_coords_sequence_error'] = '|'.join(new_sequence_errors) if any(new_sequence_errors) else ''
    allele_data.at[row_index, 'after_coords_rename_to'] = ','.join(new_allele_parts)
    allele_data.at[row_index, 'needs_fixing'] = True

allele_data.to_csv(args.output, sep='\t', index=False)

root_output_name = args.output.split('.')[0]

allele_data[allele_data['after_coords_rename_to'] != ''].to_csv(f'{root_output_name}_concerned.tsv', sep='\t', index=False)
allele_data[allele_data['after_coords_rename_to'] != ''][['allele_description', 'rename_to', 'after_coords_rename_to', 'sequence_error', 'after_coords_sequence_error']].to_csv(f'{root_output_name}_concerned_summarised.tsv', sep='\t', index=False)
