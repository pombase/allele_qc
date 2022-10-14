# %%
from genome_functions import gene_coords2genome_coords, get_feature_location_from_string
from pandas import read_csv
import warnings
from Bio import BiopythonParserWarning
import re
import pickle

data = read_csv('data/gene-coordinate-change-data.tsv', delimiter='\t', na_filter=False)


with open('data/genome.pickle', 'rb') as ins:
    contig_genome = pickle.load(ins)


def fix_coordinate(coord: str, strand):
    new_coord = coord.replace(" ", "")
    if '-' in new_coord:
        new_coord = new_coord.replace('-', '..')
    if new_coord[-1] == ',':
        new_coord = new_coord[:-1]
    if ',' not in new_coord:
        return new_coord

    match = re.search(r'\(?(\d+\.\.\d+)((,)(\d+\.\.\d+))+\)?', new_coord)

    new_coord = 'join(' + match.group().replace('(', '').replace(')', '') + ')'

    if strand == -1:
        new_coord = 'complement(' + new_coord + ')'

    return new_coord


def coord_needs_change(coord, strand):
    with warnings.catch_warnings(record=True) as w:
        try:
            loc = get_feature_location_from_string(coord)

        except AssertionError:
            print(f'error in old location in line {i + 2}: {old_coord}')

        w = list(filter(lambda i: issubclass(i.category, BiopythonParserWarning), w))
        if len(w) or (strand == -1 and 'complement' not in coord):
            return fix_coordinate(coord, strand)
        return False


for i, (systematic_id, old_coord, new_coord) in enumerate(zip(data['Systematic ID'], data['Old coordinates'], data['New coordinates'])):
    # Some of the lines have these fields empty
    if old_coord and new_coord:
        print(systematic_id)
        try:
            old_loc = get_feature_location_from_string(old_coord)
        except AssertionError:
            print(f'error in old location in line {i + 2}: {old_coord}')

        try:
            new_loc = get_feature_location_from_string(new_coord)
        except AssertionError:
            print(f'error in new location in line {i + 2}: {new_coord}')

        _, strand = gene_coords2genome_coords(1, contig_genome[systematic_id])

        c_fix = coord_needs_change(new_coord, strand)
        if c_fix:
            print(i, 'new', new_coord, c_fix, sep='\t')

        c_fix = coord_needs_change(new_coord, strand)
        if c_fix:
            print(i, 'old', old_coord, c_fix, sep='\t')
