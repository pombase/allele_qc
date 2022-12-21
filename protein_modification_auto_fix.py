import json
import pickle
import pandas
from allele_fixes import multi_shift_fix, old_coords_fix, extract_groups_from_targets, shift_coordinates_by_x


with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)

data = pandas.read_csv('results/protein_modification_results_errors_aggregated.tsv', sep='\t', na_filter=False)

with open('results/coordinate_changes_dict2.json') as ins:
    coordinate_changes_dict = json.load(ins)

# These are histone proteins that typically did not count the methionine
histones = ['SPBC1105.11c', 'SPBC1105.12', 'SPAC1834.03c', 'SPAC1834.04', 'SPAC19G12.06c', 'SPBC8D2.03c', 'SPBC8D2.04', 'SPCC622.08c', 'SPCC622.09']


def apply_multi_shift_fix(row):
    # We use at least 4 for the multi-shift
    if (row['sequence_position'].count(',') < 3):
        return ''
    if row['systematic_id'] not in genome:
        return ''
    peptide_seq = genome[row['systematic_id']]['peptide']
    return '|'.join(multi_shift_fix(peptide_seq, row['sequence_position'].split(',')))


def apply_old_coords_fix(row):
    if row['systematic_id'] not in coordinate_changes_dict:
        return '', '', ''
    result = old_coords_fix(coordinate_changes_dict[row['systematic_id']], row['sequence_position'].split(','))
    valid_solutions = pandas.DataFrame([r for r in result if all(r['matches'])])
    if valid_solutions.empty:
        return '', '', ''
    valid_solutions.loc[:, 'values'] = valid_solutions['values'].apply(','.join)
    # There could be several solutions
    return tuple('|'.join(valid_solutions[key]) for key in ['values', 'revision', 'location'])


def apply_histone_fix(row):
    if not (row['systematic_id'] in histones):
        return ''

    groups = extract_groups_from_targets(row['sequence_position'].split(','))
    peptide_seq = genome[row['systematic_id']]['peptide']

    new_coords, shifted_coords_match = shift_coordinates_by_x(groups, peptide_seq, 'number', 1)
    if all(shifted_coords_match):
        return ','.join(new_coords)

    return ''


def get_preferred_fix(row):
    """
    Based on a hierarchy
    """

    if row['histone_fix']:
        return row['histone_fix'], 'histone_fix'
    if row['old_coords_fix']:
        return row['old_coords_fix'], f'old_coords_fix, revision {row["old_coords_revision"]}: {row["old_coords_location"]}'
    if row['multi_shift_fix']:
        return row['multi_shift_fix'], 'multi_shift_fix'

    return '', ''


def format_auto_fix(row):

    syntax_error = row['change_sequence_position_to'] != ''
    sequence_position = row['sequence_position'] if not syntax_error else row['change_sequence_position_to']

    # Cases where there was no error in the sequence, but there might be a syntax error that needs fixing
    if not row['auto_fix_to']:
        return (row['change_sequence_position_to'], 'syntax_error') if syntax_error else ('', '')

    # There are edge cases where multiple fixes could be possible
    possible_fixes = set()
    for all_change_to in row['auto_fix_to'].split('|'):
        this_fix = sequence_position
        for change_from, change_to in zip(row['auto_fix_from'].split(','), all_change_to.split(',')):
            this_fix = this_fix.replace(change_from, change_to)
        possible_fixes.add(this_fix)

    return '|'.join(possible_fixes), row['auto_fix_comment']


print('applying fixes...')
extra_cols = data.apply(apply_old_coords_fix, axis=1, result_type='expand')
data.loc[:, 'old_coords_fix'] = extra_cols.iloc[:, 0]
data.loc[:, 'old_coords_revision'] = extra_cols.iloc[:, 1]
data.loc[:, 'old_coords_location'] = extra_cols.iloc[:, 2]

data['multi_shift_fix'] = data.apply(apply_multi_shift_fix, axis=1)
data['histone_fix'] = data.apply(apply_histone_fix, axis=1)

extra_cols = data.apply(get_preferred_fix, axis=1, result_type='expand')
data.loc[:, 'auto_fix_to'] = extra_cols.iloc[:, 0]
data.loc[:, 'auto_fix_comment'] = extra_cols.iloc[:, 1]
data.rename(columns={'sequence_position': 'auto_fix_from'}, inplace=True)

# Store all possible fixes
data.to_csv('results/protein_modification_auto_fix_info.tsv', sep='\t', index=False)

# Apply the fixes in the data
error_data = pandas.read_csv('results/protein_modification_results_errors.tsv', sep='\t', na_filter=False)

autofix_data = error_data.merge(data[['systematic_id', 'reference', 'auto_fix_from', 'auto_fix_to', 'auto_fix_comment']], on=['systematic_id', 'reference'], how='left')
autofix_data.fillna('', inplace=True)
extra_cols = autofix_data.apply(format_auto_fix, axis=1, result_type='expand')

# Overwrite these columns with the new values
autofix_data.loc[:, 'change_sequence_position_to'] = extra_cols.iloc[:, 0]
autofix_data.loc[:, 'auto_fix_comment'] = extra_cols.iloc[:, 1]
autofix_data.drop(columns=['auto_fix_from', 'auto_fix_to'], inplace=True)

# Print some stats
nb_errors = autofix_data.shape[0]
errors_fixed = sum(autofix_data.change_sequence_position_to != '')
syntax_errors = sum(autofix_data.auto_fix_comment == 'syntax_error')
multiple_fixes = sum(autofix_data.change_sequence_position_to.str.contains('\|'))
sequence_errors = errors_fixed - syntax_errors

print(f'{nb_errors} errors found, of which {errors_fixed} fixed:\n  - {sequence_errors} sequence errors\n  - {syntax_errors} syntax errors\n  - {multiple_fixes} have several possible fixes, for those check the `change_sequence_position_to` field for "|" characters')
print('', 'Types of errors fixed:', '', autofix_data['auto_fix_comment'].apply(lambda x: x.split(',')[0]).value_counts(), sep='\n')
# If you want to print only the dubious cases
# print(autofix_data[autofix_data.change_sequence_position_to.str.contains('\|')])
autofix_data.to_csv('results/protein_modification_auto_fix.tsv', sep='\t', index=False)
