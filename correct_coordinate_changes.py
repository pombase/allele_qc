# %%
import pandas
import re
import json

# Load alleles with errors
allele_data = pandas.read_csv('results/allele_errors.tsv', delimiter='\t', na_filter=False)
alleles_with_sequence_errors = allele_data[allele_data['sequence_error'] != '']
ids_sequence_errors = set(alleles_with_sequence_errors['systematic_id'])

# Load coordinate changes
coordinate_data = pandas.read_csv('data/all_coordinate_changes_file.tsv', delimiter='\t', na_filter=False)

# We only consider CDS features in genes that have alleles with sequence errors
coordinate_data = coordinate_data[(coordinate_data['feature_type'] == 'CDS') & (coordinate_data['systematic_id'].isin(ids_sequence_errors))]

# See the columns that only differ in value and added_removed > coordinates were modified
d = coordinate_data.drop(columns=['value', 'added_or_removed'])
logi = d.duplicated(keep=False)

# We remove here nup189, which had multiple changes
coordinate_modifications = coordinate_data[logi & (coordinate_data['primary_name'] != 'nup189')]

coordinate_modifications.to_csv('results/dummy2.tsv', sep='\t')

changes_dict = dict()

for index, row in coordinate_modifications.iterrows():
    sys_id = row['systematic_id']
    if sys_id not in changes_dict:
        changes_dict[sys_id] = dict()

    if row['added_or_removed'] == 'added':
        changes_dict[sys_id]['new'] = row['value']
    else:
        changes_dict[sys_id]['old'] = row['value']


def get_coordinate_transformer(new_coords, old_coords):
    is_complement = 'complement' in new_coords

    # THey should both have complement or not
    if is_complement != ('complement' in old_coords):
        return 'major change', ''

    # Find all numbers
    new_coords = re.findall('\d+', new_coords)
    old_coords = re.findall('\d+', old_coords)

    if len(new_coords) != len(old_coords):
        return 'major change', ''

    # Invert if they are in the -1 strand
    if is_complement:
        new_coords = new_coords[::-1]
        old_coords = old_coords[::-1]

    # Only the start has changed
    if (new_coords[0] != old_coords[0]) and (new_coords[1:] == old_coords[1:]):
        if is_complement:
            shift = int(new_coords[0]) - int(old_coords[0])
        else:
            shift = int(old_coords[0]) - int(new_coords[0])

        return 'start changed', shift
    elif (new_coords[-1] != old_coords[-1]) and (new_coords[:-1] == old_coords[:-1]):
        return 'silent change'

    return 'splicing change', ''


output_dict = dict()
not_automated = list()
for sys_id in changes_dict:
    r = changes_dict[sys_id]
    outcome, shift = get_coordinate_transformer(r['new'], r['old'])
    if shift != '':
        if shift % 3 != 0:
            print(f'error in {sys_id} shift not multiple of 3')
        output_dict[sys_id] = (0, shift)
    else:
        print(sys_id, outcome)
        print('old', r['old'])
        print('new', r['new'])

with open('results/automated_coord_fix.json', 'w') as out:
    json.dump(output_dict, out, indent=4)

with open('results/coordinate_changes.json', 'w') as out:
    json.dump(changes_dict, out, indent=4)
