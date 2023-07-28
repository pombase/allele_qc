from models import SyntaxRule, AllowedTypes
from grammar import aminoacid_grammar, allowed_types_dict, composed_types_dict, nucleotide_grammar, aminoacid_grammar_new, transition_aminoacid_grammar, transition_nucleotide_grammar, nucleotide_grammar_new
from refinement_functions import check_allele_description

import unittest

allowed_types = AllowedTypes(allowed_types=allowed_types_dict, composed_types=composed_types_dict)


class SyntaxRulesTest(unittest.TestCase):

    def test_syntax_rules(self):
        file_grammar_pairs = [
            ('test_data/aminoacid_alleles_fixable.tsv', aminoacid_grammar),
            ('test_data/aminoacid_alleles_transition.tsv', transition_aminoacid_grammar),
            ('test_data/aminoacid_alleles_new.tsv', aminoacid_grammar_new),
            ('test_data/nucleotide_alleles_fixable.tsv', nucleotide_grammar),
            ('test_data/nucleotide_alleles_transition.tsv', transition_nucleotide_grammar),
            ('test_data/nucleotide_alleles_new.tsv', nucleotide_grammar_new)
        ]
        for f, grammar in file_grammar_pairs:
            syntax_rules = [SyntaxRule.parse_obj(r) for r in grammar]

            # Remove the sequence control for testing purposes
            for r in syntax_rules:
                r.check_sequence = lambda g, gg: ''
            with open(f) as ins:
                ins.readline()
                for line_nb, line in enumerate(ins):
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

                    output = check_allele_description(allele_description, syntax_rules, allele_type, allowed_types, None)
                    try:
                        self.assertEqual(output['change_description_to'], change_description_to)
                        self.assertEqual(output['change_type_to'], change_type_to)
                        self.assertEqual(output['invalid_error'], invalid_error)
                        self.assertEqual(output['pattern_error'], pattern_error)
                    except AssertionError:
                        print(output)
                        red='\033[0;31m'
                        no_color='\033[0m'
                        print(red + f'> error in file {f} line {line_nb + 2}:' + line.strip() + no_color)
                        raise
