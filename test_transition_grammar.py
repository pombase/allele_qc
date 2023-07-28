from models import SyntaxRule
from grammar import transition_aminoacid_grammar, allowed_types, transition_nucleotide_grammar
from refinement_functions import check_allele_description

import unittest

syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in transition_aminoacid_grammar]
syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in transition_nucleotide_grammar]

# Remove the sequence control for testing purposes
for r in syntax_rules_aminoacids:
    r.check_sequence = lambda g, gg: ''

# Remove the sequence control for testing purposes
for r in syntax_rules_nucleotides:
    r.check_sequence = lambda g, gg: ''


class TransitionGrammarTest(unittest.TestCase):

    def test_aminoacid_cases(self):
        with open('test_data/aminoacid_alleles_transition.tsv') as ins:
            ins.readline()
            for line in ins:
                if len(line.strip()) == 0 or line[0] == '#':
                    continue
                ls = line.strip().split('\t')
                allele_type, allele_description, change_description_to, change_type_to, invalid_error = [ls[0], ls[1], '', '', '']
                if len(ls) > 2:
                    change_description_to = ls[2]
                if len(ls) > 3:
                    change_type_to = ls[3]
                if len(ls) > 4:
                    invalid_error = ls[4]
                output = check_allele_description(allele_description, syntax_rules_aminoacids, allele_type, allowed_types, None)
                try:
                    self.assertEqual(output['change_description_to'], change_description_to)
                    self.assertEqual(output['change_type_to'], change_type_to)
                    self.assertEqual(output['invalid_error'], invalid_error)
                except AssertionError:
                    print('error in line:', line.strip())
                    raise

    def test_nucleotide_cases(self):
        with open('test_data/nucleotide_alleles_transition.tsv') as ins:
            ins.readline()
            for line in ins:
                if len(line.strip()) == 0 or line[0] == '#':
                    continue
                ls = line.strip().split('\t')
                allele_type, allele_description, change_description_to, change_type_to, invalid_error, pattern_error = [ls[0], ls[1], '', '', '', '']
                if len(ls) > 2:
                    change_description_to = ls[2]
                if len(ls) > 3:
                    change_type_to = ls[3]
                if len(ls) > 4:
                    invalid_error = ls[4]
                if len(ls) > 5:
                    pattern_error = ls[5]

                output = check_allele_description(allele_description, syntax_rules_nucleotides, allele_type, allowed_types, None)
                try:
                    self.assertEqual(output['change_description_to'], change_description_to)
                    self.assertEqual(output['change_type_to'], change_type_to)
                    self.assertEqual(output['invalid_error'], invalid_error)
                    self.assertEqual(output['pattern_error'], pattern_error)
                except AssertionError:
                    print('error in line:', line.strip())
                    raise
