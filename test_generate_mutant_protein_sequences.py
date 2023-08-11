import unittest
from generate_mutant_protein_sequences import transvar_variant_to_substitution_dict, variant_sequence_from_subsitution_dicts


class GenerateMutantProteinTest(unittest.TestCase):

    def test_generate_mutant_sequence(self):

        sequence = 'MAVLQQQQ*'

        # Test single substitution
        variant = 'p.A2B'
        expected = {'range': slice(1, 2), 'replace_by': 'B'}
        self.assertEqual(transvar_variant_to_substitution_dict(variant, len(sequence)), expected)
        expected_sequence = 'MBVLQQQQ*'
        self.assertEqual(variant_sequence_from_subsitution_dicts(sequence, [expected]), expected_sequence)

        # Test partial deletion
        variant = 'p.A2_Q5delAVLQ'
        expected = {'range': slice(1, 5), 'replace_by': ''}
        self.assertEqual(transvar_variant_to_substitution_dict(variant, len(sequence)), expected)
        expected_sequence = 'MQQQ*'
        self.assertEqual(variant_sequence_from_subsitution_dicts(sequence, [expected]), expected_sequence)

        # Test single deletion
        variant = 'p.A2delA'
        expected = {'range': slice(1, 2), 'replace_by': ''}
        self.assertEqual(transvar_variant_to_substitution_dict(variant, len(sequence)), expected)
        expected_sequence = 'MVLQQQQ*'
        self.assertEqual(variant_sequence_from_subsitution_dicts(sequence, [expected]), expected_sequence)

        # Test insertion
        variant = 'p.A2_V3insCCCCC'
        expected = {'range': slice(2, 2), 'replace_by': 'CCCCC'}
        self.assertEqual(transvar_variant_to_substitution_dict(variant, len(sequence)), expected)
        expected_sequence = 'MACCCCCVLQQQQ*'
        self.assertEqual(variant_sequence_from_subsitution_dicts(sequence, [expected]), expected_sequence)

        # Should give an error if the start and stop are not adjacent
        variant = 'p.A2_L4insCCCCC'
        try:
            transvar_variant_to_substitution_dict(variant, len(sequence))
        except ValueError:
            pass
        else:
            raise AssertionError('Should have raised an error')

        # Test duplication
        variant = 'p.A2_L4dupAVL'
        expected = {'range': slice(1, 1), 'replace_by': 'AVL'}
        self.assertEqual(transvar_variant_to_substitution_dict(variant, len(sequence)), expected)
        expected_sequence = 'MAVLAVLQQQQ*'
        self.assertEqual(variant_sequence_from_subsitution_dicts(sequence, [expected]), expected_sequence)

        # Should give an error if the duplicated sequence does not match the length
        variant = 'p.A2_Q5dupAVL'
        try:
            transvar_variant_to_substitution_dict(variant, len(sequence))
        except ValueError:
            pass
        else:
            raise AssertionError('Should have raised an error')

        # Test delins
        variant = 'p.A2_Q5delinsCCC'
        expected = {'range': slice(1, 5), 'replace_by': 'CCC'}
        self.assertEqual(transvar_variant_to_substitution_dict(variant, len(sequence)), expected)
        expected_sequence = 'MCCCQQQ*'
        self.assertEqual(variant_sequence_from_subsitution_dicts(sequence, [expected]), expected_sequence)

        # Test multiple substitutions
        variants = ['p.A2B', 'p.V3C', 'p.V3_L4insCCCCC', 'p.Q5_Q6delinsTTTT', 'p.Q7_Q8delQQ']
        variant_list = [transvar_variant_to_substitution_dict(variant, len(sequence)) for variant in variants]
        expected = [{'range': slice(1, 2), 'replace_by': 'B'},
                    {'range': slice(2, 3), 'replace_by': 'C'},
                    {'range': slice(3, 3), 'replace_by': 'CCCCC'},
                    {'range': slice(4, 6), 'replace_by': 'TTTT'},
                    {'range': slice(6, 8), 'replace_by': ''}]
        self.assertEqual(variant_list, expected)
        sequence = 'MAVLQQQQ*'
        expected_sequence = 'MBCCCCCCLTTTT*'
        self.assertEqual(variant_sequence_from_subsitution_dicts(sequence, variant_list), expected_sequence)

        # Test multiple insertions passed in either order
        variants = ['p.A2_V3insCC', 'p.V3_L4insTT']
        variant_list = [transvar_variant_to_substitution_dict(variant, len(sequence)) for variant in variants]
        expected = [{'range': slice(2, 2), 'replace_by': 'CC'},
                    {'range': slice(3, 3), 'replace_by': 'TT'}]
        self.assertEqual(variant_list, expected)
        sequence = 'MAVLQQQQ*'
        expected_sequence = 'MACCVTTLQQQQ*'
        self.assertEqual(variant_sequence_from_subsitution_dicts(sequence, variant_list), expected_sequence)
        self.assertEqual(variant_sequence_from_subsitution_dicts(sequence, variant_list[::-1]), expected_sequence)

        # Test stop codon
        variant = 'p.L4*'
        sequence = 'MAVLQQQQ*'
        expected = {'range': slice(3, 9), 'replace_by': '*'}
        self.assertEqual(transvar_variant_to_substitution_dict(variant, len(sequence)), expected)
        expected_sequence = 'MAV*'
        self.assertEqual(variant_sequence_from_subsitution_dicts(sequence, [expected]), expected_sequence)
