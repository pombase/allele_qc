from Bio import pairwise2
import re
from genome_functions import get_other_index_from_alignment

new_seq = 'MNKSIFIQLQDQIDKEHSIREKLTAEVDLLDEKLRVLQLLLANCEQNLENQEEILEALEIIKSKTRGLAELASNFPYYKYNGVWDRSIQKVVYLYLLASWTGRLDKSLRPTYSLLSLSEVGQILQVPVFPEESTFHLSIEQYLHAVLSLCSELARQSVNSVISGNYHIPFEALNTIQKVHSSFQVLSLKNDSLRRHFDGLKYDLKRSEDVVYDLRIHKLV*'
old_seq = 'MNKSIFIQLQDQIDKEHSIREKLTAEVDLLDEKLRVLQLLLANCEQSRNENLQEKEHGLTLEDLENQEEILEALEIIKSKTRGLAELASNFPYYKYNGVWDRSIQKVVYLYLLASWTGRLDKSLRPTYSLLSLSEVGQILQVPVFPEESTFHLSIEQYLHAVLSLCSELARQSVNSVISGNYHIPFEALNTIQKVHSSFQVLSLKNDSLRRHFDGLKYDLKRSEDVVYDLRIHKLV*'
# An alignment where gaps have the maximum penalty, to avoid scattered matches.
alignments = pairwise2.align.globalms(new_seq, old_seq, match=2, mismatch=-1, open=-100, extend=0, penalize_end_gaps=False)

new_alignment = alignments[0].seqA
old_alignment = alignments[0].seqB

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
