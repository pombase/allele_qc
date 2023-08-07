from models import SyntaxRule, AllowedTypes, find_rule
from grammar import aminoacid_grammar_old, allowed_types_dict, composed_types_dict, nucleotide_grammar_old,\
    aminoacid_grammar, transition_old2new_aminoacid_grammar, transition_old2new_nucleotide_grammar, nucleotide_grammar,\
    transition_new2old_aminoacid_grammar, transition_new2old_nucleotide_grammar

from refinement_functions import check_allele_description
import unittest
import pickle


allowed_types = AllowedTypes(allowed_types=allowed_types_dict, composed_types=composed_types_dict)

with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)


class SyntaxRulesTest(unittest.TestCase):

    def test_class_methods_amino_acids(self):
        syntax_rules = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]

        syntax_rule = find_rule(syntax_rules, 'amino_acid_mutation', 'multiple_aa')
        self.assertEqual(syntax_rule.type, 'amino_acid_mutation')
        self.assertEqual(syntax_rule.rule_name, 'multiple_aa')
        groups = syntax_rule.get_groups('AP123VL')
        self.assertEqual(groups, ('AP', '123', 'VL'))
        self.assertEqual(syntax_rule.format_for_transvar(groups), 'p.A123_P124delinsVL')

        syntax_rule = find_rule(syntax_rules, 'amino_acid_insertion', 'standard')
        self.assertEqual(syntax_rule.type, 'amino_acid_insertion')
        self.assertEqual(syntax_rule.rule_name, 'standard')
        groups = syntax_rule.get_groups('E2EVTA')
        self.assertEqual(groups, ('E', '2', 'EVTA'))
        self.assertEqual(syntax_rule.format_for_transvar(groups), 'p.E2_3insEVTA')

        syntax_rule = find_rule(syntax_rules, 'nonsense_mutation', 'stop_codon_star')
        self.assertEqual(syntax_rule.type, 'nonsense_mutation')
        self.assertEqual(syntax_rule.rule_name, 'stop_codon_star')
        groups = syntax_rule.get_groups('E2*')
        self.assertEqual(groups, ('E', '2', '*'))
        self.assertEqual(syntax_rule.format_for_transvar(groups), 'p.E2*')

        syntax_rule = find_rule(syntax_rules, 'partial_amino_acid_deletion', 'multiple_aa')
        self.assertEqual(syntax_rule.type, 'partial_amino_acid_deletion')
        self.assertEqual(syntax_rule.rule_name, 'multiple_aa')
        groups = syntax_rule.get_groups('2-100')
        self.assertEqual(groups, ('2', '100'))
        self.assertEqual(syntax_rule.format_for_transvar(groups), 'p.2_100del')

        syntax_rule = find_rule(syntax_rules, 'partial_amino_acid_deletion', 'single_aa')
        self.assertEqual(syntax_rule.type, 'partial_amino_acid_deletion')
        self.assertEqual(syntax_rule.rule_name, 'single_aa')
        groups = syntax_rule.get_groups('2')
        self.assertEqual(groups, ('2',))
        self.assertEqual(syntax_rule.format_for_transvar(groups), 'p.2del')

    def test_class_methods_nucleotides(self):
        syntax_rules = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]

        # Gene in the +1 strand
        ase1_gene = genome['SPAPB1A10.09']
        syntax_rule = find_rule(syntax_rules, 'nucleotide_mutation', 'multiple_nt')
        self.assertEqual(syntax_rule.type, 'nucleotide_mutation')
        self.assertEqual(syntax_rule.rule_name, 'multiple_nt')
        # Should be Q2M at the protein level
        groups = syntax_rule.get_groups('CAA4ATG')
        self.assertEqual(groups, ('CAA', '4', 'ATG'))
        self.assertEqual(syntax_rule.format_for_transvar(groups, ase1_gene), 'g.1878365_1878367delCAAinsATG')

        # In the 5UTR
        groups = syntax_rule.get_groups('CAA(-4)ATG')
        self.assertEqual(groups, ('CAA', '(-4)', 'ATG'))
        # self.assertEqual(syntax_rule.format_for_transvar(groups, ase1_gene), 'g.1878365_1878367delCAAinsATG')

        # Gene in the -1 strand
        mse1_gene = genome['SPAPB1A10.11c']
        syntax_rule = find_rule(syntax_rules, 'nucleotide_mutation', 'multiple_nt')
        self.assertEqual(syntax_rule.type, 'nucleotide_mutation')
        self.assertEqual(syntax_rule.rule_name, 'multiple_nt')
        # Should be T5Y at the protein level
        groups = syntax_rule.get_groups('ACC4TAT')
        self.assertEqual(groups, ('ACC', '4', 'TAT'))
        self.assertEqual(syntax_rule.format_for_transvar(groups, mse1_gene), 'g.1884046_1884048delGGTinsATA')

    def test_syntax_rules(self):
        file_grammar_pairs = [
            ('test_data/aminoacid_alleles_fixable.tsv', aminoacid_grammar_old),
            ('test_data/aminoacid_alleles_transition_old2new.tsv', transition_old2new_aminoacid_grammar),
            ('test_data/aminoacid_alleles_transition_new2old.tsv', transition_new2old_aminoacid_grammar),
            ('test_data/aminoacid_alleles_new.tsv', aminoacid_grammar),
            ('test_data/nucleotide_alleles_fixable.tsv', nucleotide_grammar_old),
            ('test_data/nucleotide_alleles_transition_old2new.tsv', transition_old2new_nucleotide_grammar),
            ('test_data/nucleotide_alleles_transition_new2old.tsv', transition_new2old_nucleotide_grammar),
            ('test_data/nucleotide_alleles_new.tsv', nucleotide_grammar)
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
