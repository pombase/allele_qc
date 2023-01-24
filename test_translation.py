"""
Check that sequence translations match the ones in pombase
"""
import unittest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import requests
import gzip
import io
import pickle

with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)


known_exceptions = ['SPAC977.01']


class TranslationTest(unittest.TestCase):

    def test_translations(self):
        # Get peptide data from pombase
        response = requests.get('https://www.pombase.org/releases/latest/fasta/feature_sequences/peptide.fa.gz')
        data = gzip.decompress(response.content).decode()
        for seq in SeqIO.parse(io.StringIO(data), 'fasta'):
            seq: SeqRecord

            # A peptide or an intron
            if seq.id.count('.') > 1:
                systematic_id = '.'.join(seq.id.split('.')[:2])
            else:
                systematic_id = seq.id

            if systematic_id not in genome or 'peptide' not in genome[systematic_id] or systematic_id in known_exceptions:
                continue

            self.assertEqual(genome[systematic_id]['peptide'], seq.seq)
