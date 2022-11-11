from Bio import pairwise2
import re
from genome_functions import get_other_index_from_alignment, get_feature_location_from_string
import pickle

new_seq = 'MNKSIFIQLQDQIDKEHSIREKLTAEVDLLDEKLRVLQLLLANCEQNLENQEEILEALEIIKSKTRGLAELASNFPYYKYNGVWDRSIQKVVYLYLLASWTGRLDKSLRPTYSLLSLSEVGQILQVPVFPEESTFHLSIEQYLHAVLSLCSELARQSVNSVISGNYHIPFEALNTIQKVHSSFQVLSLKNDSLRRHFDGLKYDLKRSEDVVYDLRIHKLV*'
old_seq = 'MNKSIFIQLQDQIDKEHSIREKLTAEVDLLDEKLRVLQLLLANCEQSRNENLQEKEHGLTLEDLENQEEILEALEIIKSKTRGLAELASNFPYYKYNGVWDRSIQKVVYLYLLASWTGRLDKSLRPTYSLLSLSEVGQILQVPVFPEESTFHLSIEQYLHAVLSLCSELARQSVNSVISGNYHIPFEALNTIQKVHSSFQVLSLKNDSLRRHFDGLKYDLKRSEDVVYDLRIHKLV*'
# An alignment where gaps have the maximum penalty, to avoid scattered matches.
alignments = pairwise2.align.globalms(new_seq, old_seq, match=2, mismatch=-1, open=-100, extend=0, penalize_end_gaps=False)

new_alignment = alignments[0].seqA
old_alignment = alignments[0].seqB

print(old_alignment)
print(new_alignment)


targets = ['S204stop',
           'Y227stop',
           'K106A',
           'K95A',
           'Q105A',
           'R102G',
           'R20G',
           'R210G',
           'R211G',
           'R35G',
           'W100A', ]

numbers = list()
values = list()
for t in targets:
    numbers.append(int(re.search(r'\d+', t).group()) - 1)
    values.append(t[0])

for v, old_index in zip(values, numbers):
    new_index = get_other_index_from_alignment(old_alignment, new_alignment, old_index)

    print(v, old_seq[old_index], old_index + 1, new_seq[new_index], new_index + 1, sep='\t')

# %%
with open('data/genome.pickle', 'rb') as ins:
    contig_genome = pickle.load(ins)

new_feature_loc = get_feature_location_from_string('complement(join(492370..492784,492832..493062,493108..493549,493600..493738,493780..493954,494003..495612,495673..495722,495767..495926))')
old_feature_loc = get_feature_location_from_string('complement(join(492370..492784,492832..493549,493600..493738,493780..493954,494003..495612,495673..495722,495767..495926))')

new_seq = new_feature_loc.extract(contig_genome['SPCC970.09']['contig']).translate()
old_seq = old_feature_loc.extract(contig_genome['SPCC970.09']['contig']).translate()

# An alignment where gaps have the maximum penalty, to avoid scattered matches.
alignments = pairwise2.align.globalms(new_seq.seq, old_seq.seq, match=2, mismatch=-1, open=-100, extend=0, penalize_end_gaps=False)

new_alignment = alignments[0].seqA
old_alignment = alignments[0].seqB


targets = ['Q992*', ]

numbers = list()
values = list()
for t in targets:
    numbers.append(int(re.search(r'\d+', t).group()) - 1)
    values.append(t[0])

for v, old_index in zip(values, numbers):
    new_index = get_other_index_from_alignment(old_alignment, new_alignment, old_index)

    print(v, old_seq[old_index], old_index + 1, new_seq[new_index], new_index + 1, sep='\t')
