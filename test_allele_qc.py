import unittest
from allele_qc import handle_systematic_id_for_qc
import pandas

class HandleSystematicIdTest(unittest.TestCase):

    def test_handle_systematic_id(self):

        genome = {
            'SPBC1198.04c.1': {'CDS': 'dummy', 'peptide': 'MSAS'},
            'SPBC1198.04c.2': {'CDS': 'dummy', 'peptide': 'AAAAAA'},
            'dummy': {'CDS': 'dummy', 'peptide': 'MSAS'}
        }

        row = pandas.Series({'systematic_id': 'SPBC1198.04c', 'gene_name': 'zas1', 'allele_name': 'zas1.2-A1V'})
        self.assertEqual(handle_systematic_id_for_qc(row, genome), 'SPBC1198.04c.2')

        row = pandas.Series({'systematic_id': 'SPBC1198.04c', 'gene_name': 'zas1', 'allele_name': 'zas1.1-A1V'})
        self.assertEqual(handle_systematic_id_for_qc(row, genome), 'SPBC1198.04c.1')

        row = pandas.Series({'systematic_id': 'SPBC1198.04c', 'gene_name': 'zas1', 'allele_name': 'zas1-A1V'})
        self.assertEqual(handle_systematic_id_for_qc(row, genome), 'SPBC1198.04c.1')

        row = pandas.Series({'systematic_id': 'dummy', 'gene_name': 'blah', 'allele_name': 'blah-A1V'})
        self.assertEqual(handle_systematic_id_for_qc(row, genome), 'dummy')


