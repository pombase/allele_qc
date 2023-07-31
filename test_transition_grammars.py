from models import SyntaxRule, AllowedTypes
from grammar import aminoacid_grammar, allowed_types_dict, composed_types_dict, nucleotide_grammar, transition_old2new_aminoacid_grammar, transition_old2new_nucleotide_grammar, disruption_grammar
import pandas
from allele_qc import check_fun
import unittest
import pickle


class TransitionGrammarsTest(unittest.TestCase):
    # To test that conversion is reversible

    def test_transition_grammar(self):
        allowed_types = AllowedTypes(allowed_types=allowed_types_dict, composed_types=composed_types_dict)

        with open('data/genome.pickle', 'rb') as ins:
            genome = pickle.load(ins)
        allele_data = pandas.read_csv('data/alleles.tsv', delimiter='\t', na_filter=False)
        syntax_rules_aminoacids_old2new = [SyntaxRule.parse_obj(r) for r in transition_old2new_aminoacid_grammar]
        syntax_rules_nucleotides_old2new = [SyntaxRule.parse_obj(r) for r in transition_old2new_nucleotide_grammar]
        syntax_rules_disruption = [SyntaxRule.parse_obj(r) for r in disruption_grammar]
        allowed_types = AllowedTypes(allowed_types=allowed_types_dict, composed_types=composed_types_dict)

        extra_cols = allele_data.apply(lambda row: check_fun(row, genome, syntax_rules_aminoacids_old2new, syntax_rules_nucleotides_old2new, syntax_rules_disruption, allowed_types), axis=1, result_type='expand')
        new_fixes = pandas.concat([allele_data, extra_cols], axis=1)
        # Keep only those with corrections
        new_fixes = new_fixes[(new_fixes['change_description_to'] != '') & (new_fixes['pattern_error'] == '') & (new_fixes['invalid_error'] == '')]
        print(new_fixes)
