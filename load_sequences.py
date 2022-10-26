"""
Build a genome dictionary from several fasta files. It might be PomBase-specific. We favour the
use of these sequences for peptide sequences since they take into account nucleotide changes in
the genome sequence that have not been incorporated to the contig files. The dictionary is stored
as a pickle file, as in 'load_genome.py'

There is a fasta file for each type of sequences:

data
├── cds+introns+utrs.fa
├── cds+introns.fa
├── cds.fa
├── five_prime_utrs.fa
├── introns_within_cds.fa
├── peptide.fa
└── three_prime_utrs.fa

And each file contains the sequences of that kind, with the systematic_id, e.g. for peptide.

>SPAC1002.01.1:pep mrx11|mitochondrial expression network (MIOREX) component Mrx11
MLPPTIRISGLAKTLHIPSRSPLQALKGSFILLNKRKFHYSPFILQEKVQSSNHTIRSDT
KLWKRLLKITGKQAHQFKDKPFSHIFAFLFLHELSAILPLPIFFFIFHSLDWTPTGLPGE
YLQKGSHVAASIFAKLGYNLPLEKVSKTLLDGAAAYAVVKVC*
>SPAC1002.02.1:pep pom34|nucleoporin Pom34
MASTFSQSVFARSLYEDSAENKVDSSKNTEANFPITLPKVLPTDPKASSLHKPQEQQPNI
IPSKEEDKKPVINSMKLPSIPAPGTDNINESHIPRGYWKHPAVDKIAKRLHDQAPSDRTW
SRMVSNLFAFISIQFLNRYLPNTTAVKVVSWILQALLLFNLLESVWQFVRPQPTFDDLQL
TPLQRKLMGLPEGGSTSGKHLTPPRYRPNFSPSRKAENVKSPVRSTTWA*

The script generates a dictionary in which the keys are the systematic_id of genes and the values are dictionaries
where the keys are the titles of the fasta files, and the value is the sequence. E.g.:

{
    "SPAC1002.01.1": {
        "peptide": "MLPPTIRISG...",
        "CDS": "ATG...",
        etc.
    }
}
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import pickle


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('files', metavar='N', type=str, nargs='+',
                    help='files to be read')
parser.add_argument('--output', default='data/fasta_genome.pickle')
args = parser.parse_args()

fasta_genome = dict()

for f in args.files:
    print('\033[0;32m reading: ' + f + '\033[0m')
    feature_type = f[5:].split('.')[0]
    for seq in SeqIO.parse(f, 'fasta'):
        seq: SeqRecord

        # A peptide or an intron
        if seq.id.count('.') > 1:
            systematic_id = '.'.join(seq.id.split('.')[:2])
        else:
            systematic_id = seq.id

        if systematic_id not in fasta_genome:
            fasta_genome[systematic_id] = dict()

        fasta_genome[systematic_id][feature_type] = seq.seq

with open(args.output, 'wb') as out:
    pickle.dump(fasta_genome, out, pickle.HIGHEST_PROTOCOL)
