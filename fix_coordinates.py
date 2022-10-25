from models import SyntaxRule, find_rule
import argparse
import pandas
import pickle
from genome_functions import get_other_index_from_alignment
from refinement_functions import get_allele_parts_from_result, build_regex2syntax_rule, replace_allele_features, sort_result
import re
import json
from grammar import aminoacid_grammar

parser = argparse.ArgumentParser(description='Build a dictionary of alignments based on the updated coordinates of genes')
parser.add_argument('--genome', default='data/genome.pickle')
parser.add_argument('--alleles', default='results/allele_errors.tsv')
parser.add_argument('--coordinate_changes', default='results/coordinate_changes_dict.json')
parser.add_argument('--alleles_excluded', default='results/ambiguous_coordinate_changes.tsv')
args = parser.parse_args()


with open(args.genome, 'rb') as ins:
    contig_genome = pickle.load(ins)

with open(args.coordinate_changes, 'r') as ins:
    changes_dict = json.load(ins)

# Load alleles with sequence errors
allele_data = pandas.read_csv(args.alleles, delimiter='\t', na_filter=False)
allele_data = allele_data[(allele_data['sequence_error'] != '') & allele_data['systematic_id'].isin(changes_dict) & ~allele_data['allele_type'].str.contains('nucl')]

# Exclude ambiguous
alleles_excluded = pandas.read_csv(args.alleles_excluded, delimiter='\t', na_filter=False)
allele_data = allele_data[~allele_data['allele_name'].isin(alleles_excluded['allele_name'])]

syntax_rules = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
regex2syntax_rule = build_regex2syntax_rule(syntax_rules)

for i, row in allele_data.iterrows():

    rules_in_this_row = [find_rule(syntax_rules, rule_type, rule_name) for rule_type, rule_name in [rule_pair.split(':') for rule_pair in row['rules_applied'].split('|')]]
    new_alignment = changes_dict[row['systematic_id']]['new_alignment']
    old_alignment = changes_dict[row['systematic_id']]['old_alignment']

    new_allele_parts = list()
    for rule, allele_part in zip(rules_in_this_row, row['allele_parts'].split('|')):
        match = re.match(rule.regex, allele_part)
        groups = match.groups()
        for i in rule.coordinate_indexes:
            new_index = get_other_index_from_alignment(int(groups[i]))
            groups[i] = str(new_index)

        if any(g == 'None' for g in groups):
            break
        rule.check_sequence(groups, fasta_genome[systematic_id])

    # regex2syntax_rule = build_regex2syntax_rule(syntax_rules)

    # result = replace_allele_features(list(regex2syntax_rule.keys()), [row['allele_description']], [])

    # allele_parts = get_allele_parts_from_result(result)

    # # Extract the matched and unmatched elements
    # matches, unmatched = sort_result(result)

    # new_alignment = changes_dict[row['systematic_id']]['new_alignment']
    # old_alignment = changes_dict[row['systematic_id']]['old_alignment']

    # for i, match in enumerate(matches):
    #     syntax_rule = regex2syntax_rule[match.re.pattern]
    #     old_groups = match.groups()
    #     new_allele_descritpion
