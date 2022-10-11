# %%
from genome_functions import get_feature_location_from_string
from pandas import read_csv


data = read_csv('data/gene-coordinate-change-data.tsv', delimiter='\t', na_filter=False)

for i, (systematic_id, old_coord, new_coord) in enumerate(zip(data['Systematic ID'], data['Old coordinates'], data['New coordinates'])):
    # Some of the lines have these fields empty
    if old_coord and new_coord:
        try:
            old_loc = get_feature_location_from_string(old_coord)
        except AssertionError:
            print(f'error in old location in line {i + 2}: {old_coord}')

        try:
            new_loc = get_feature_location_from_string(new_coord)
        except AssertionError:
            print(f'error in new location in line {i + 2}: {new_coord}')

        # new_loc = get_feature_location_from_string(new_coord)

        # except AssertionError:
        #     print(old_coord)
        #     print(new_coord)
