from grammar import grammar
from models import SyntaxRule
from grammar import grammar, allowed_types
from refinement_functions import find_allele_parts

import unittest

syntax_rules = [SyntaxRule.parse_obj(r) for r in grammar]

# Remove the sequence control for testing purposes
for r in syntax_rules:
    r.check_sequence = lambda g, gg: ''


class SyntaxRulesTest(unittest.TestCase):

    def test_first(self):
        with open('test_data/alleles_fixable.tsv') as ins:
            ins.readline()
            for line in ins:
                if len(line.strip()) == 0 or line[0] == '#':
                    continue
                ls = line.strip().split('\t')
                allele_type, allele_description, rename_to, change_type_to, invalid_error = [ls[0], ls[1], '', '', '']
                if len(ls) > 2:
                    rename_to = ls[2]
                if len(ls) > 3:
                    change_type_to = ls[3]
                if len(ls) > 4:
                    invalid_error = ls[4]
                output = find_allele_parts(allele_description, syntax_rules, allele_type, allowed_types, None)
                self.assertEqual(output['rename_to'], rename_to)
                self.assertEqual(output['change_type_to'], change_type_to)
                self.assertEqual(output['invalid_error'], invalid_error)
