import unittest
from genome_functions import get_other_index_from_alignment
from allele_fixes import multi_shift_fix, old_coords_fix


class SequenceIndexingTest(unittest.TestCase):

    def test_coordinate_changes(self):

        examples = (
            (
                "VAQCAAIKVT----AQCVKVTVIFLAAAA",
                "VAQC--IKVTVIFLAQCVKVTVIFL"
            ),
            (
                "VAQC--IKVTVIFLAQCVKVTVIFL",
                "VAQCAAIKVT----AQCVKVTVIFLAAAA"
            ),
            (
                "-VAT---",
                "AVATCGF"
            ),
            (
                "AVATCGF",
                "-VAT---"
            )
        )

        for new_alignment, old_alignment in examples:
            new_seq = new_alignment.replace('-', '')
            old_seq = old_alignment.replace('-', '')
            for i in range(len(new_seq)):
                old_i = get_other_index_from_alignment(new_alignment, old_alignment, i)
                if old_i is None:
                    continue
                self.assertEqual(new_seq[i], old_seq[old_i])

    def test_shift_coordinates(self):
        seq = 'AVPPAVPPP'
        self.assertEqual(multi_shift_fix(seq, ['A3', 'V4']), ['A1,V2', 'A5,V6'])
        self.assertEqual(multi_shift_fix(seq, ['A3', 'V4', '8-10']), ['A1,V2,6-8'])

        # Test that it finds both extremes
        self.assertEqual(multi_shift_fix('AVVVVVVVVA', ['A3']), ['A1', 'A10'])

        # Test that it works with one-length string
        self.assertEqual(multi_shift_fix('V', ['A3']), [])
        self.assertEqual(multi_shift_fix('A', ['A3']), ['A1'])

    def test_old_coords_fix(self):

        coordinate_changes = [
            # Correct
            {
                'revision': 'dummy_revision',
                'new_alignment': 'ACLPTV',
                'old_alignment': 'A--PTV',
                'old_coord': 'dummy_coord'
            },
            # Old does not match
            {
                'revision': 'dummy_revision2',
                'new_alignment': 'ACLPTV',
                'old_alignment': 'V--PTV',
                'old_coord': 'dummy_coord1'
            },
            # Another correct one (to see that it returns both)
            {
                'revision': 'dummy_revision3',
                'new_alignment': 'ACPTV',
                'old_alignment': 'A-PTV',
                'old_coord': 'dummy_coord3'
            },
            # New does not match
            {
                'revision': 'dummy_revision4',
                'new_alignment': 'VCLPTV',
                'old_alignment': 'A--PTV',
                'old_coord': 'dummy_coord4'
            },
            # Position does not exist in old
            {
                'revision': 'dummy_revision5',
                'new_alignment': 'VCLPTV',
                'old_alignment': 'A-----',
                'old_coord': 'dummy_coord5'
            },
            # Position does not exist in new
            {
                'revision': 'dummy_revision6',
                'new_alignment': '-CLPTV',
                'old_alignment': 'A--PTV',
                'old_coord': 'dummy_coord6'
            },
        ]

        solutions = old_coords_fix(coordinate_changes, ['A1', 'P2'])
        self.assertEqual(len(solutions), 2)
        self.assertEqual(solutions[0]['values'], 'A1,P4')
        self.assertEqual(solutions[1]['values'], 'A1,P3')
