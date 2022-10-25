import unittest
from genome_functions import get_other_index_from_alignment


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
