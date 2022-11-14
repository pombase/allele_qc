from Bio import pairwise2
from genome_functions import get_other_index_from_alignment, get_feature_location_from_string
import pickle
import pandas

modif_coords = pandas.read_csv('data/only_modified_coordinates.tsv', sep='\t', na_filter=False)

systematic_id = 'SPAC630.14c'
modif_coords = modif_coords[modif_coords.systematic_id == systematic_id]

# The last coordinate is the last one with 'added'
last_coord = list(modif_coords['value'])[0]

prev_coords = list(modif_coords['value'][1:][modif_coords.added_or_removed == 'removed'])

with open('data/genome.pickle', 'rb') as ins:
    contig_genome = pickle.load(ins)

for prev_coord in prev_coords:

    new_feature_loc = get_feature_location_from_string(last_coord)
    old_feature_loc = get_feature_location_from_string(prev_coord)

    new_seq = new_feature_loc.extract(contig_genome[systematic_id]['contig']).translate()
    old_seq = old_feature_loc.extract(contig_genome[systematic_id]['contig']).translate()

    # An alignment where gaps have the maximum penalty, to avoid scattered matches.
    alignments = pairwise2.align.globalms(new_seq.seq, old_seq.seq, match=2, mismatch=-1, open=-100, extend=0, penalize_end_gaps=False)

    new_alignment = alignments[0].seqA
    old_alignment = alignments[0].seqB

    new_indexes = list()
    old_indexes = [67, 271]
    for old_index in old_indexes:
        new_indexes.append(get_other_index_from_alignment(old_alignment, new_alignment, old_index))

    print(len(old_seq), ': ', *old_indexes, '->', *new_indexes, prev_coord)
