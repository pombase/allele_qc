"""
Build a dictionary of alignments based on the updated coordinates of genes, and store it as json.

The input is a coordinate changes file generated with https://github.com/pombase/genome_changelog

The dictionary structure, where keys are the systematic_id of genes:

"SPAC23E2.02": [{
        "new_revision": "new_revision",
        "old_revision": "old_revision",
        "new_coord": "join(446770..449241,449295..450530)",
        "old_coord": "join(446491..446513,446679..449241,449295..450530)",
        "new_alignment": "--------------------------------------MNTSENDP ... GYNGTRY*",
        "old_alignment": "MPLGRSSWICCAKYFVNTKSRFNEILPPRFTLIVSFYSMNTSENDP ... SGYNGTRY*"
}],

In the alignment gaps have the maximum penalty, to avoid scattered matches.

-----CAV
MAACATAV

And not

---C--AV
MAACATAV

The reason for this is that the alignment is based on changing / removing introns or changing the start of ending
coordinate of the start or end of the CDS, so you want maximal identity with minimum number of gaps.
"""

import argparse
import pandas
import pickle
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from genome_functions import get_feature_location_from_string
import json
from pathlib import Path
import glob


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--genome', default='data/genome.pickle')
parser.add_argument('--coords', default='data/only_modified_coordinates.tsv')
parser.add_argument('--output', default='results/coordinate_changes_dict2.json')
args = parser.parse_args()

chromosome_dictionary = {
    'I': 'chromosome1',
    'II': 'chromosome2',
    'III': 'chromosome3',
    'mating_type_region': 'mating_type_region',
    'mitochondrial': 'pMIT',
}


def read_old_genomes(files, format):
    """
    Create a dictionary of dictionaries where the keys are contig_name and revision and value is
    the sequence of the given contig at a given revision, such that to access the sequence of
    chromosome X at revision Y, you do output[x][y].

    It relies on a file structure where files are like this:
    ├── chromosome1
    │   ├── 20011004.contig
    │   ├── 20020322.contig
    ....
    """
    out_dict = dict()
    for f in files:
        splitted_path = Path(f).parts
        revision = splitted_path[-1].split('.')[-2]
        contig = splitted_path[-2]

        if contig not in out_dict:
            out_dict[contig] = dict()

        with open(f, errors='replace') as ins:
            # We store it as seq because we don't want to store
            out_dict[contig][revision] = Seq(SeqIO.read(ins, format)[0])

    return out_dict


def choose_old_genome(previous_coordinate, latest_genome_seq, old_genomes_dict, genome_seq_changes: pandas.DataFrame):
    if any(genome_seq_changes[['chromosome', 'date']].duplicated()):
        raise ValueError('The script cannot handle two revisions made on the same date')

    # Sort by date
    genome_seq_changes = genome_seq_changes.sort_values(['date'], ascending=False).copy()

    # Select which genome should be used
    contig = chromosome_dictionary[previous_coordinate['chromosome']]
    changes_this_contig = genome_seq_changes[genome_seq_changes['chromosome'] == contig]
    newer_genome_sequences = previous_coordinate['date'] <= changes_this_contig['date']

    # All revisions where changes were made are newer, use the oldest genome (they are sorted)
    chromosome_revisions_dictionary = old_genomes_dict[contig]
    if all(newer_genome_sequences):
        return SeqRecord(chromosome_revisions_dictionary[changes_this_contig['previous_revision'].iloc[-1]])
    # All revisions where changes were made are older, use the newest genome
    elif all(~newer_genome_sequences):
        return latest_genome_seq

    # previous_revision of the next version
    chosen_revision = changes_this_contig.loc[newer_genome_sequences, 'previous_revision'].iloc[-1]
    return SeqRecord(chromosome_revisions_dictionary[chosen_revision])


# Load info about changes in genome sequence
genome_seq_changes = pandas.read_csv('data/genome_sequence_changes.tsv', sep='\t', na_filter=False, dtype=str)

print('reading old genomes...')
old_genomes_dict = read_old_genomes(glob.glob('data/old_genome_versions/*/*.contig'), 'embl')
print('old genomes read')

with open(args.genome, 'rb') as ins:
    latest_genome = pickle.load(ins)

# Load coordinate changes
coordinate_data = pandas.read_csv(args.coords, delimiter='\t', na_filter=False)

# We only consider CDS features in genes that have alleles with sequence errors
coordinate_data = coordinate_data[(coordinate_data['feature_type'] == 'CDS')].copy()

# Make sure the order is right
coordinate_data.sort_values(['date', 'revision', 'chromosome', 'systematic_id', 'feature_type', 'added_or_removed'], inplace=True, ascending=[False, False, True, True, True, True])

changes_dict = dict()

for systematic_id in set(coordinate_data.systematic_id):
    data_subset = coordinate_data[coordinate_data.systematic_id == systematic_id].copy()
    # First added is the latest coordinates
    latest_coordinate = data_subset[data_subset['added_or_removed'] == 'added'].iloc[0, :]
    # The rest of "removed" are the older ones
    previous_coordinates = data_subset[data_subset['added_or_removed'] == 'removed']
    changes_dict[systematic_id] = list()

    for i, previous_coordinate in previous_coordinates.iterrows():
        this_change = dict()
        this_change['revision'] = previous_coordinate['revision']
        this_change['new_coord'] = latest_coordinate['value']
        this_change['old_coord'] = previous_coordinate['value']

        new_feature_loc = get_feature_location_from_string(latest_coordinate['value'])
        old_feature_loc = get_feature_location_from_string(previous_coordinate['value'])
        if systematic_id not in latest_genome or 'contig' not in latest_genome[systematic_id]:
            print(systematic_id, 'skipped, probably to do with alt-splicing')
            continue

        new_seq = new_feature_loc.extract(latest_genome[systematic_id]['contig']).translate()
        # The genome on which a feature was defined might not have the same sequence, if so use the old sequence
        old_genome = choose_old_genome(previous_coordinate, latest_genome[systematic_id]['contig'], old_genomes_dict, genome_seq_changes)
        old_seq = old_feature_loc.extract(old_genome).translate()

        # This can happen when the sequence was changed more than once to the same thing.
        if new_seq.seq == old_seq.seq:
            continue

        # An alignment where gaps have the maximum penalty, to avoid scattered matches.
        alignments = pairwise2.align.globalms(new_seq.seq, old_seq.seq, match=1, mismatch=-2, open=-2, extend=0, penalize_end_gaps=False)
        if len(alignments) == 0:
            continue

        this_change['new_alignment'] = alignments[0].seqA
        this_change['old_alignment'] = alignments[0].seqB

        changes_dict[systematic_id].append(this_change)

with open(args.output, 'w') as out:
    json.dump(changes_dict, out, indent=4)
