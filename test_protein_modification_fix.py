from common_autofix_functions import format_auto_fix
import unittest

from protein_modification_qc import check_func


class ProteinModificationFormattingTest(unittest.TestCase):

    def test_correct_syntax_error(self):
        genome = {
            'dummy': {'CDS': 'dummy', 'peptide': 'MSAS'}
        }
        example_dict = {
            'systematic_id': 'dummy',
            'sequence_position': 'S4',
        }
        self.assertEqual(check_func(example_dict, genome), ('', ''))

        example_dict['sequence_position'] = 'S3'
        self.assertEqual(check_func(example_dict, genome), ('S3', ''))

        example_dict['sequence_position'] = 'blah3'
        self.assertEqual(check_func(example_dict, genome), ('pattern_error', ''))

        example_dict['sequence_position'] = 'A3S4'
        self.assertEqual(check_func(example_dict, genome), ('pattern_error', ''))

        example_dict['sequence_position'] = 'A3,S4'
        self.assertEqual(check_func(example_dict, genome), ('', ''))

        # Repetitions are kept (expected behaviour)
        example_dict['sequence_position'] = 'A3,S4,S4'
        self.assertEqual(check_func(example_dict, genome), ('', ''))

        # Only wrong syntax errors are reported, '|' separated
        example_dict['sequence_position'] = 'A3,S4,T4,T1000'
        self.assertEqual(check_func(example_dict, genome), ('||T4|T1000', ''))

        example_dict['sequence_position'] = 'A3;S4'
        self.assertEqual(check_func(example_dict, genome), ('', 'A3,S4'))

        # No sorting, but repetitions are kept even when correcting syntax error
        example_dict['sequence_position'] = 'S4 A3 S4;S4;T100'
        self.assertEqual(check_func(example_dict, genome), ('||||T100', 'S4,A3,S4,S4,T100'))

    def test_format_auto_fix(self):

        # Simplest case
        example_dict = {
            'sequence_position': 'S409',
            'sequence_error': 'S409',
            'change_sequence_position_to': '',
            'auto_fix_from': 'S409',
            'auto_fix_to': 'S400',
            'auto_fix_comment': 'blah'
        }
        fixes, comment = format_auto_fix(example_dict, 'sequence_position', 'change_sequence_position_to')
        self.assertEqual(fixes, 'S400')
        self.assertEqual(comment, 'blah')

        # There was an error sometime in which I was using string.replace
        # This gave an error as it would pick S4 from S409 and give S109
        # (replacing the substring S4 by S1)
        example_dict['auto_fix_from'] = 'S4,S409'
        example_dict['auto_fix_to'] = 'S1,S400'
        fixes, comment = format_auto_fix(example_dict, 'sequence_position', 'change_sequence_position_to')
        self.assertEqual(fixes, 'S400')
        self.assertEqual(comment, 'blah')

        # Case where both old solutions make no difference
        example_dict['auto_fix_to'] = 'S1,S400|S3,S400'
        example_dict['auto_fix_comment'] = 'blah|bluh'
        fixes, comment = format_auto_fix(example_dict, 'sequence_position', 'change_sequence_position_to')
        self.assertEqual(fixes, 'S400')
        self.assertEqual(comment, 'blah/bluh')

        # Case where both solutions make a difference
        example_dict['auto_fix_to'] = 'S1,S401|S3,S400'
        fixes, comment = format_auto_fix(example_dict, 'sequence_position', 'change_sequence_position_to')
        self.assertEqual(fixes, 'S401|S400')
        self.assertEqual(comment, 'blah|bluh')

        # Autofix works with random strings as well
        example_dict['auto_fix_to'] = 'S1,??'
        fixes, comment = format_auto_fix(example_dict, 'sequence_position', 'change_sequence_position_to')
        self.assertEqual(fixes, '??')
