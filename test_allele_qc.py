import unittest
from allele_qc import handle_systematic_id_for_allele_qc
import pickle


with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)


class HandleSystematicIdTest(unittest.TestCase):

    def test_handle_systematic_id(self):

        row = {'systematic_id': 'SPBC1198.04c', 'allele_name': 'zas1.2-A1V'}
        self.assertEqual(handle_systematic_id_for_allele_qc(row['systematic_id'], row['allele_name'], genome), 'SPBC1198.04c.2')

        row = {'systematic_id': 'SPBC1198.04c', 'allele_name': 'zas1.1-A1V'}
        self.assertEqual(handle_systematic_id_for_allele_qc(row['systematic_id'], row['allele_name'], genome), 'SPBC1198.04c.1')

        row = {'systematic_id': 'SPBC1198.04c', 'allele_name': 'zas1-A1V'}
        self.assertEqual(handle_systematic_id_for_allele_qc(row['systematic_id'], row['allele_name'], genome), 'SPBC1198.04c.1')



