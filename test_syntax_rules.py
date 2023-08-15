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

    def test_find_syntax_rule(self):
        syntax_rules = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
        syntax_rule = find_rule(syntax_rules, 'amino_acid_mutation', 'multiple_aa')
        self.assertEqual(syntax_rule.type, 'amino_acid_mutation')
        self.assertEqual(syntax_rule.rule_name, 'multiple_aa')

        syntax_rules = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]
        syntax_rule = find_rule(syntax_rules, 'nucleotide_mutation', 'multiple_nt')
        self.assertEqual(syntax_rule.type, 'nucleotide_mutation')
        self.assertEqual(syntax_rule.rule_name, 'multiple_nt')

    def test_class_methods_amino_acids(self):
        syntax_rules = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]

        syntax_rule = find_rule(syntax_rules, 'amino_acid_mutation', 'multiple_aa')
        groups = syntax_rule.get_groups('AP123VL', None)
        self.assertEqual(groups, ('AP', '123', 'VL'))
        self.assertEqual(syntax_rule.format_for_transvar(groups, None), 'p.A123_P124delinsVL')

        syntax_rule = find_rule(syntax_rules, 'amino_acid_insertion', 'standard')
        groups = syntax_rule.get_groups('E2EVTA', None)
        self.assertEqual(groups, ('E', '2', 'EVTA'))
        self.assertEqual(syntax_rule.format_for_transvar(groups, None), 'p.E2_3insVTA')

        syntax_rule = find_rule(syntax_rules, 'nonsense_mutation', 'stop_codon_star')
        groups = syntax_rule.get_groups('E2*', None)
        self.assertEqual(groups, ('E', '2', '*'))
        self.assertEqual(syntax_rule.format_for_transvar(groups, None), 'p.E2*')

        syntax_rule = find_rule(syntax_rules, 'partial_amino_acid_deletion', 'multiple_aa')
        groups = syntax_rule.get_groups('2-100', None)
        self.assertEqual(groups, ('2', '100'))
        self.assertEqual(syntax_rule.format_for_transvar(groups, None), 'p.2_100del')

        syntax_rule = find_rule(syntax_rules, 'partial_amino_acid_deletion', 'single_aa')
        groups = syntax_rule.get_groups('2', None)
        self.assertEqual(groups, ('2',))
        self.assertEqual(syntax_rule.format_for_transvar(groups, None), 'p.2del')

    def test_class_methods_CTD(self):
        syntax_rules = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]

        syntax_rule = find_rule(syntax_rules, 'amino_acid_mutation', 'CTD')
        groups = syntax_rule.get_groups('CTD-Y1F(r1-r12),Y2F', genome['SPBC28F2.12'])
        self.assertEqual(groups, ('CTD-Y1F(r1-r12),Y2F', ))

        syntax_rule = find_rule(syntax_rules, 'partial_amino_acid_deletion', 'CTD')
        groups = syntax_rule.get_groups('CTD-delta', genome['SPBC28F2.12'])
        self.assertEqual(groups, ('CTD-delta', ))

        syntax_rule = find_rule(syntax_rules, 'partial_amino_acid_deletion', 'CTD')
        groups = syntax_rule.get_groups('CTD-delta(r1-r12-2)', genome['SPBC28F2.12'])
        self.assertEqual(groups, ('CTD-delta(r1-r12-2)', ))

        syntax_rule = find_rule(syntax_rules, 'partial_amino_acid_deletion', 'CTD')
        groups = syntax_rule.get_groups('CTD-Δ(r1-r12-2)', genome['SPBC28F2.12'])
        self.assertEqual(groups, ('CTD-Δ(r1-r12-2)', ))

        syntax_rule = find_rule(syntax_rules, 'amino_acid_deletion_and_mutation', 'CTD')
        groups = syntax_rule.get_groups('CTD-delta(r1-r12-2),Y1F(r1-r12),Y2F', genome['SPBC28F2.12'])
        self.assertEqual(groups, ('CTD-delta(r1-r12-2),Y1F(r1-r12),Y2F', ))

    def test_class_methods_nucleotides(self):

        syntax_rules = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]
        syntax_rule_multi_nt = find_rule(syntax_rules, 'nucleotide_mutation', 'multiple_nt')

        # Gene in the +1 strand
        ase1_gene = genome['SPAPB1A10.09']
        # Should be Q2M at the protein level
        groups = syntax_rule_multi_nt.get_groups('CAA4ATG', None)
        self.assertEqual(groups, ('CAA', '4', 'ATG'))
        self.assertEqual(syntax_rule_multi_nt.format_for_transvar(groups, ase1_gene), 'g.1878365_1878367delCAAinsATG')

        # In the 5UTR (negative number)
        groups = syntax_rule_multi_nt.get_groups('GTTCA(-5)CCCAC', None)
        self.assertEqual(groups, ('GTTCA', '(-5)', 'CCCAC'))
        self.assertEqual(syntax_rule_multi_nt.format_for_transvar(groups, ase1_gene), 'g.1878357_1878361delGTTCAinsCCCAC')

        # Gene in the -1 strand
        mse1_gene = genome['SPAPB1A10.11c']
        # Should be T5Y at the protein level
        groups = syntax_rule_multi_nt.get_groups('ACC4TAT', None)
        self.assertEqual(groups, ('ACC', '4', 'TAT'))
        self.assertEqual(syntax_rule_multi_nt.format_for_transvar(groups, mse1_gene), 'g.1884046_1884048delGGTinsATA')

        # In the 5UTR (negative number)
        groups = syntax_rule_multi_nt.get_groups('TTTTGG(-6)CCCAAA', None)
        self.assertEqual(groups, ('TTTTGG', '(-6)', 'CCCAAA'))
        self.assertEqual(syntax_rule_multi_nt.format_for_transvar(groups, mse1_gene), 'g.1884052_1884057delCCAAAAinsTTTGGG')

        # Single nucleotide substitution
        syntax_rule_single_nt = find_rule(syntax_rules, 'nucleotide_mutation', 'single_nt')

        groups = syntax_rule_single_nt.get_groups('A1G', None)
        self.assertEqual(groups, ('A', '1', 'G'))
        self.assertEqual(syntax_rule_single_nt.format_for_transvar(groups, ase1_gene), 'g.1878362A>G')

        groups = syntax_rule_single_nt.get_groups('A(-1)G', None)
        self.assertEqual(groups, ('A', '(-1)', 'G'))
        self.assertEqual(syntax_rule_single_nt.format_for_transvar(groups, ase1_gene), 'g.1878361A>G')

        groups = syntax_rule_single_nt.get_groups('A1G', None)
        self.assertEqual(groups, ('A', '1', 'G'))
        self.assertEqual(syntax_rule_single_nt.format_for_transvar(groups, mse1_gene), 'g.1884051T>C')

        groups = syntax_rule_single_nt.get_groups('G(-1)A', None)
        self.assertEqual(groups, ('G', '(-1)', 'A'))
        self.assertEqual(syntax_rule_single_nt.format_for_transvar(groups, mse1_gene), 'g.1884052C>T')

        # Nucleotide insertion
        syntax_rule_insertion = find_rule(syntax_rules, 'nucleotide_insertion', 'standard')

        groups = syntax_rule_insertion.get_groups('A1AGGG', None)
        self.assertEqual(groups, ('A', '1', 'AGGG'))
        self.assertEqual(syntax_rule_insertion.format_for_transvar(groups, ase1_gene), 'g.1878362_1878363insGGG')

        groups = syntax_rule_insertion.get_groups('A(-1)AGGG', None)
        self.assertEqual(groups, ('A', '(-1)', 'AGGG'))
        self.assertEqual(syntax_rule_insertion.format_for_transvar(groups, ase1_gene), 'g.1878361_1878362insGGG')

        groups = syntax_rule_insertion.get_groups('A1AGGG', None)
        self.assertEqual(groups, ('A', '1', 'AGGG'))
        self.assertEqual(syntax_rule_insertion.format_for_transvar(groups, mse1_gene), 'g.1884050_1884051insCCC')

        groups = syntax_rule_insertion.get_groups('G(-1)GTTT', None)
        self.assertEqual(groups, ('G', '(-1)', 'GTTT'))
        self.assertEqual(syntax_rule_insertion.format_for_transvar(groups, mse1_gene), 'g.1884051_1884052insAAA')

        # Single nucleotide deletion
        syntax_rule_deletion = find_rule(syntax_rules, 'partial_nucleotide_deletion', 'single_nt')

        groups = syntax_rule_deletion.get_groups('1', None)
        self.assertEqual(groups, ('1',))
        self.assertEqual(syntax_rule_deletion.format_for_transvar(groups, ase1_gene), 'g.1878362del')

        groups = syntax_rule_deletion.get_groups('-1', None)
        self.assertEqual(groups, ('-1',))
        self.assertEqual(syntax_rule_deletion.format_for_transvar(groups, ase1_gene), 'g.1878361del')

        # Multi nucleotide deletion
        syntax_rule_deletion = find_rule(syntax_rules, 'partial_nucleotide_deletion', 'usual')

        groups = syntax_rule_deletion.get_groups('1-10', None)
        self.assertEqual(groups, ('1', '10'))
        self.assertEqual(syntax_rule_deletion.format_for_transvar(groups, ase1_gene), 'g.1878362_1878371del')

        # Note that this function would not work if the order was inverted, so we should only run on descriptions
        # that have been corrected if needed.
        groups = syntax_rule_deletion.get_groups('(-10)-(-1)', None)
        self.assertEqual(groups, ('(-10)', '(-1)'))
        self.assertEqual(syntax_rule_deletion.format_for_transvar(groups, ase1_gene), 'g.1878352_1878361del')

        groups = syntax_rule_deletion.get_groups('1-10', None)
        self.assertEqual(groups, ('1', '10'))
        self.assertEqual(syntax_rule_deletion.format_for_transvar(groups, mse1_gene), 'g.1884042_1884051del')

        groups = syntax_rule_deletion.get_groups('(-10)-(-1)', None)
        self.assertEqual(groups, ('(-10)', '(-1)'))
        self.assertEqual(syntax_rule_deletion.format_for_transvar(groups, mse1_gene), 'g.1884052_1884061del')

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

                    # We use genome['SPBC28F2.12'] for it to work with the CTD ones
                    output = check_allele_description(allele_description, syntax_rules, allele_type, allowed_types, genome['SPBC28F2.12'])
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
