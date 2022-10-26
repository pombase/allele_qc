"""
These tests only work with PomBase data.
"""
import unittest
import pickle
from grammar import get_nt_at_genome_position
from Bio.SeqIO import parse

with open('data/genome.pickle', 'rb') as ins:
    contig_genome = pickle.load(ins)

with open('data/fasta_genome.pickle', 'rb') as ins:
    fasta_genome = pickle.load(ins)


class SequenceIndexingTest(unittest.TestCase):

    def test_nts(self):
        # Coding and non-coding examples, from the +1 and -1 strands
        test_genes = [
            {
                'id': 'SPAPB1A10.09',
                'downstream': 'ATGCAAACAGTAATGATGG',
                'upstream': 'TCATTTACATCAACCGGTTCA',
            },
            {
                'id': 'SPAPB1A10.10c',
                'downstream': 'ATGTCGGCTCAGAAAAGGG',
                'upstream': 'AGGTACGACAGAATATACTTCA',
            },
            {
                'id': 'SPNCRNA.2846',
                'downstream': 'ACTTCTTTTTGCTTGCAAAGTT',
                'upstream': 'CTTTTCTTTTTCAGCGGAAAAA',
            },
            {
                'id': 'SPNCRNA.2847',
                'downstream': 'ACATTATAAACAATTACACAACAATCGGCCCCTC',
                'upstream': 'ACTGAGTCAAAAGACTTCGAGTTATTC',
            },
        ]
        for g in test_genes:
            # Gene from the forward strand (ase1)
            gene = contig_genome[g['id']]

            # Check that downstream matches
            for i, value in enumerate(g['downstream']):
                self.assertEqual(get_nt_at_genome_position(i + 1, gene, gene['contig']), value)

            # Check that 5'UTR matches
            for i, value in enumerate(g['upstream'][::-1]):
                self.assertEqual(get_nt_at_genome_position(-i, gene, gene['contig']), value)

    def test_aas(self):
        for seq in parse('test_data/test_peptides.fasta', 'fasta'):
            self.assertEqual(fasta_genome[seq.id]['peptide'], seq.seq)
            self.assertEqual(contig_genome[seq.id]['peptide'], seq.seq)
