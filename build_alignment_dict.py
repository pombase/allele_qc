import argparse
import pandas
import pickle
from Bio import pairwise2
from genome_functions import get_feature_location_from_string
import json

parser = argparse.ArgumentParser(description='Build a dictionary of alignments based on the updated coordinates of genes')
parser.add_argument('--genome', default='data/genome.pickle')
parser.add_argument('--alleles', default='results/allele_errors.tsv')
parser.add_argument('--coords', default='data/all_coordinate_changes_file.tsv')
parser.add_argument('--output', default='results/coordinate_changes_dict.json')
args = parser.parse_args()


with open(args.genome, 'rb') as ins:
    contig_genome = pickle.load(ins)


# Load alleles with errors
allele_data = pandas.read_csv(args.alleles, delimiter='\t', na_filter=False)
alleles_with_sequence_errors = allele_data[allele_data['sequence_error'] != '']
ids_sequence_errors = set(alleles_with_sequence_errors['systematic_id'])

# Load coordinate changes
coordinate_data = pandas.read_csv(args.coords, delimiter='\t', na_filter=False)

# We only consider CDS features in genes that have alleles with sequence errors
coordinate_data = coordinate_data[(coordinate_data['feature_type'] == 'CDS') & (coordinate_data['systematic_id'].isin(ids_sequence_errors))]

# See the columns that only differ in value and added_removed > coordinates were modified
d = coordinate_data.drop(columns=['value', 'added_or_removed'])
logi = d.duplicated(keep=False)

# We remove here nup189, which had multiple changes
coordinate_modifications = coordinate_data[logi & (coordinate_data['primary_name'] != 'nup189')]

changes_dict = dict()

for index, row in coordinate_modifications.iterrows():
    sys_id = row['systematic_id']
    if sys_id not in changes_dict:
        changes_dict[sys_id] = dict()

    if row['added_or_removed'] == 'added':
        changes_dict[sys_id]['new_coord'] = row['value']
    else:
        changes_dict[sys_id]['old_coord'] = row['value']

for sys_id in changes_dict:
    new_coord = changes_dict[sys_id]['new_coord']
    old_coord = changes_dict[sys_id]['old_coord']

    new_feature_loc = get_feature_location_from_string(new_coord)
    old_feature_loc = get_feature_location_from_string(old_coord)
    if 'contig' not in contig_genome[sys_id]:
        print(sys_id, 'skipped, probably to do with alt-splicing')
        continue

    new_seq = new_feature_loc.extract(contig_genome[sys_id]['contig']).translate()
    old_seq = old_feature_loc.extract(contig_genome[sys_id]['contig']).translate()

    # An alignment where gaps have the maximum penalty, to avoid scattered matches.
    alignments = pairwise2.align.globalms(new_seq.seq, old_seq.seq, match=2, mismatch=-1, open=-100, extend=0, penalize_end_gaps=False)

    changes_dict[sys_id]['new_alignment'] = alignments[0].seqA
    changes_dict[sys_id]['old_alignment'] = alignments[0].seqB

with open(args.output, 'w') as out:
    json.dump(changes_dict, out, indent=4)
