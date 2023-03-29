import json
import pickle
import pandas
from common_autofix_functions import apply_multi_shift_fix, apply_old_coords_fix, apply_histone_fix, get_preferred_fix, format_auto_fix


with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)

data = pandas.read_csv('results/protein_modification_results_errors_aggregated.tsv', sep='\t', na_filter=False)

with open('data/coordinate_changes_dict.json') as ins:
    coordinate_changes_dict = json.load(ins)


print('applying fixes...')
extra_cols = data.apply(apply_old_coords_fix, axis=1, result_type='expand', args=[coordinate_changes_dict, 'sequence_position'])
data.loc[:, 'old_coords_fix'] = extra_cols.iloc[:, 0]
data.loc[:, 'old_coords_revision'] = extra_cols.iloc[:, 1]
data.loc[:, 'old_coords_location'] = extra_cols.iloc[:, 2]

data['multi_shift_fix'] = data.apply(apply_multi_shift_fix, axis=1, args=[genome, 'sequence_position'])
data['histone_fix'] = data.apply(apply_histone_fix, axis=1, args=[genome, 'sequence_position'])

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
extra_cols = autofix_data.apply(format_auto_fix, axis=1, result_type='expand', args=['sequence_position', 'change_sequence_position_to'])

# Overwrite these columns with the new values
autofix_data.loc[:, 'change_sequence_position_to'] = extra_cols.iloc[:, 0].apply(lambda x: x.split('|'))
autofix_data.loc[:, 'auto_fix_comment'] = extra_cols.iloc[:, 1].apply(lambda x: x.split('|'))
autofix_data.drop(columns=['auto_fix_from', 'auto_fix_to'], inplace=True)

# Explode columns with multiple solutions
autofix_data.loc[:, 'solution_index'] = autofix_data['change_sequence_position_to'].apply(lambda x: list(range(len(x))) if len(x) > 1 else [None, ])
autofix_data = autofix_data.explode(['change_sequence_position_to', 'auto_fix_comment', 'solution_index'])

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

fixed_rows = autofix_data.change_sequence_position_to != ''
autofix_data[fixed_rows].to_csv('results/protein_modification_auto_fix.tsv', sep='\t', index=False)

other_errors_names = ['systematic_id not in genome', 'missing_CDS', 'pattern_error']

cannot_fix = autofix_data[~fixed_rows].drop(columns=['change_sequence_position_to', 'auto_fix_comment', 'solution_index'])

other_errors = cannot_fix.sequence_error.isin(other_errors_names)

cannot_fix[other_errors].rename(columns={'sequence_error': 'error'}).to_csv('results/protein_modification_cannot_fix_other_errors.tsv', sep='\t', index=False)
cannot_fix[~other_errors].to_csv('results/protein_modification_cannot_fix_sequence_errors.tsv', sep='\t', index=False)
