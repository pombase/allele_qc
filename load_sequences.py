# %%
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob

fasta_files = glob.glob('data/*.fa')

fasta_genome = dict()

for f in fasta_files:
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
