import pandas
from grammar import aminoacid_grammar
from models import SyntaxRule
from refinement_functions import split_multiple_aa, join_multiple_aa
import pickle
import json
import re
from common_autofix_functions import apply_multi_shift_fix, apply_old_coords_fix, apply_histone_fix, get_preferred_fix
import os


with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)

with open('results/coordinate_changes_dict.json') as ins:
    coordinate_changes_dict = json.load(ins)

syntax_rules = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
syntax_rules_dict = {f'{r.type}:{r.rule_name}': r for r in syntax_rules}

data = pandas.read_csv('results/allele_results_errors.tsv', sep='\t', na_filter=False)

# We only want alleles with sequence errors, with mutations described at the aminoacid level
aminoacid_alleles_with_sequence_errors = (data['sequence_error'] != '') & (data['rules_applied'].str.contains('amino_acid') | data['rules_applied'].str.contains('nonsense'))
data_subset = data.loc[aminoacid_alleles_with_sequence_errors, ['systematic_id', 'allele_description', 'allele_name', 'reference', 'change_description_to', 'rules_applied', 'sequence_error']]

# Explode the references
data_subset.loc[:, 'reference'] = data_subset['reference'].apply(str.split, args=[','])
data_subset = data_subset.explode('reference')

# Keep correct description only
data_subset.loc[data_subset['change_description_to'] != '', 'allele_description'] = data_subset.loc[data_subset['change_description_to'] != '', 'change_description_to']
data_subset.drop(columns='change_description_to', inplace=True)

# Copy allele_description to explode it
data_subset['allele_description_exploded'] = data_subset['allele_description'].copy()

# Explode
cols2explode = ['allele_description_exploded', 'rules_applied', 'sequence_error']
for explode_column, sep in zip(cols2explode, [',', '|', '|']):
    data_subset.loc[:, explode_column] = data_subset[explode_column].apply(str.split, args=[sep])
data_subset = data_subset.explode(cols2explode)


# Explode again the multiple_aa
multi_aa = data_subset['rules_applied'] == 'amino_acid_mutation:multiple_aa'
data_subset['allele_description_exploded2'] = data_subset['allele_description_exploded'].copy()
data_subset.loc[multi_aa, 'allele_description_exploded2'] = data_subset.loc[multi_aa, 'allele_description_exploded2'].apply(split_multiple_aa, args=[syntax_rules_dict['amino_acid_mutation:multiple_aa'].regex])
data_subset = data_subset.explode('allele_description_exploded2')

# Sort the mutations by the first number in them (there might be two numbers in the deletions)
data_subset.loc[:, 'sorting_col'] = data_subset['allele_description_exploded2'].apply(lambda x: int(re.search(r'\d+', x).group()))
data_subset.sort_values('sorting_col', inplace=True)
data_subset.drop(columns='sorting_col', inplace=True)

data_subset.drop_duplicates(inplace=True)


# Make a copy for matching the proposed fixes later
data_for_fixing = data_subset.copy()

# Drop duplicates that may have arised when splitting allele description at either step
data_subset.drop(columns=['rules_applied', 'allele_description', 'allele_description_exploded'], inplace=True)
data_subset.drop_duplicates(inplace=True)


# Aggregate by publication and systematic id
aggregated_data = data_subset.groupby(['systematic_id', 'reference'], as_index=False).agg({'allele_description_exploded2': ','.join})


print('applying fixes...')
extra_cols = aggregated_data.apply(apply_old_coords_fix, axis=1, result_type='expand', args=[coordinate_changes_dict, 'allele_description_exploded2'])
aggregated_data.loc[:, 'old_coords_fix'] = extra_cols.iloc[:, 0]
aggregated_data.loc[:, 'old_coords_revision'] = extra_cols.iloc[:, 1]
aggregated_data.loc[:, 'old_coords_location'] = extra_cols.iloc[:, 2]
aggregated_data['multi_shift_fix'] = aggregated_data.apply(apply_multi_shift_fix, axis=1, args=[genome, 'allele_description_exploded2'])
aggregated_data['histone_fix'] = aggregated_data.apply(apply_histone_fix, axis=1, args=[genome, 'allele_description_exploded2'])
extra_cols = aggregated_data.apply(get_preferred_fix, axis=1, result_type='expand')
aggregated_data.loc[:, 'auto_fix_to'] = extra_cols.iloc[:, 0]
aggregated_data.loc[:, 'auto_fix_comment'] = extra_cols.iloc[:, 1]
aggregated_data.to_csv('results/allele_auto_fix_info.tsv', sep='\t', index=False)

# Re-explode the columns that have multiple solutions (now they are aggregated as '|'-separated strings)
aggregated_data_with_fixes = aggregated_data.loc[aggregated_data['auto_fix_to'] != '', :].copy()

# Deaggregate multiple solutions
aggregated_data_with_fixes.loc[:, 'auto_fix_to'] = aggregated_data_with_fixes['auto_fix_to'].apply(str.split, args=['|'])

# We add an extra column with the index of the solution if there is more than one
# We convert to strings because the missing values may give unexpected behaviour in .groupby or .agg
aggregated_data_with_fixes.loc[:, 'solution_index'] = aggregated_data_with_fixes['auto_fix_to'].apply(lambda x: list(range(len(x))) if len(x) > 1 else [''])
aggregated_data_with_fixes.loc[:, 'solution_index'] = aggregated_data_with_fixes.loc[:, 'solution_index']
aggregated_data_with_fixes = aggregated_data_with_fixes.explode(['auto_fix_to', 'solution_index'])

# Deaggregate the parts of each solution
aggregated_data_with_fixes.loc[:, 'allele_description_exploded2'] = aggregated_data_with_fixes['allele_description_exploded2'].apply(str.split, args=[','])
aggregated_data_with_fixes.loc[:, 'auto_fix_to'] = aggregated_data_with_fixes['auto_fix_to'].apply(str.split, args=[','])

for i, row in aggregated_data_with_fixes.iterrows():
    if len(row['allele_description_exploded2']) != len(row['auto_fix_to']):
        print(row['allele_description_exploded2'], row['auto_fix_to'])

data2merge = aggregated_data_with_fixes.explode(['allele_description_exploded2', 'auto_fix_to'])
data_for_fixing = data_for_fixing.merge(data2merge[['systematic_id', 'reference', 'allele_description_exploded2', 'auto_fix_to', 'auto_fix_comment', 'solution_index']], on=['systematic_id', 'reference', 'allele_description_exploded2'])

# Aggregate the multiple_aa, they only differ on the columns allele_description_exploded2 and auto_fix_to, so we aggregate on everything else
# We also drop allele_description_exploded and allele_description_exploded2 and rules_applied after the aggregation, no longer needed
multi_aa = data_for_fixing['rules_applied'] == 'amino_acid_mutation:multiple_aa'
groupby_columns = list(data_for_fixing.columns)
groupby_columns.remove('allele_description_exploded2')
groupby_columns.remove('auto_fix_to')
new_rows = data_for_fixing[multi_aa].groupby(groupby_columns, as_index=False).agg({'auto_fix_to': lambda x: join_multiple_aa(x.values)})
data_for_fixing = data_for_fixing[~multi_aa]
data_for_fixing = pandas.concat([data_for_fixing, new_rows])
data_for_fixing.drop(columns=['allele_description_exploded', 'allele_description_exploded2', 'sequence_error', 'rules_applied'], inplace=True)
data_for_fixing.drop_duplicates(inplace=True)

# Sort again by first number in the string (concat step may mess up the order,
# which is important for the next aggregation)
data_for_fixing.loc[:, 'sorting_col'] = data_for_fixing['auto_fix_to'].apply(lambda x: int(re.search(r'\d+', x).group()))
data_for_fixing.sort_values('sorting_col', inplace=True)
data_for_fixing.drop(columns='sorting_col', inplace=True)

# Apply the fix by merging the auto_fix_to individual columns
groupby_columns = list(data_for_fixing.columns)
groupby_columns.remove('auto_fix_to')

data_for_fixing = data_for_fixing.groupby(groupby_columns, as_index=False).agg({'auto_fix_to': ','.join})

# Merge solutions from different PMIDs that are the same
groupby_columns = list(data_for_fixing.columns)
groupby_columns.remove('reference')
data_for_fixing = data_for_fixing.groupby(groupby_columns, as_index=False).agg({'reference': ','.join})

# Here we check if multiple solutions were found from different references, if so give a warning
groupby_columns.remove('auto_fix_to')
different_solutions = data_for_fixing[data_for_fixing[groupby_columns].duplicated(keep=False)]
if not different_solutions.empty:
    print('Different solutions have been found in different papers')
    print(different_solutions)

# Finally, we merge with the original data
data = data.merge(data_for_fixing[['systematic_id', 'allele_name', 'auto_fix_comment', 'solution_index', 'auto_fix_to']], on=['systematic_id', 'allele_name'], how='outer')
data.fillna('', inplace=True)

# Set auto_fix_to value in change_description_to, and drop the column
fixed_sequence_errors = data['auto_fix_to'] != ''
data.loc[fixed_sequence_errors, 'change_description_to'] = data.loc[fixed_sequence_errors, 'auto_fix_to']
data.drop(columns='auto_fix_to', inplace=True)

columns_auto_fixed = fixed_sequence_errors | \
    ((data['sequence_error'] == '') & ((data['change_type_to'] != '') | (data['change_description_to'] != '')))

autofixed_data = data[columns_auto_fixed].copy()

# Fill comments for the rest of error types
syntax_error = (autofixed_data.auto_fix_comment == '') & (autofixed_data.change_description_to != '')
type_error = (autofixed_data.auto_fix_comment == '') & (autofixed_data.change_type_to != '')

autofixed_data.loc[syntax_error, 'auto_fix_comment'] = 'syntax_error'
autofixed_data.loc[type_error, 'auto_fix_comment'] = 'type_error'
autofixed_data.loc[syntax_error & type_error, 'auto_fix_comment'] = 'syntax_and_type_error'

autofixed_data = autofixed_data[['systematic_id', 'allele_description', 'allele_name', 'allele_type', 'change_description_to', 'change_type_to', 'auto_fix_comment', 'sequence_error', 'solution_index', 'allele_parts', 'rules_applied', 'reference']].copy()

autofixed_data.to_csv('results/allele_auto_fix.tsv', sep='\t', index=False)

errors_cannot_fix = data[~columns_auto_fixed].drop(columns=['auto_fix_comment', 'solution_index'])

# If there is a file with manual fixes, do not include in cannot_fix
if os.path.isfile('manual_fixes_pombase/manual_changes_formatted.tsv'):
    manual_changes = pandas.read_csv('manual_fixes_pombase/manual_changes_formatted.tsv', sep='\t', na_filter=False)
    manual_changes['combined_column'] = manual_changes.apply(lambda r: [r['systematic_id'], r['allele_name']], axis=1)
    errors_cannot_fix['combined_column'] = errors_cannot_fix.apply(lambda r: [r['systematic_id'], r['allele_name']], axis=1)
    errors_cannot_fix = errors_cannot_fix[~errors_cannot_fix['combined_column'].isin(manual_changes['combined_column'])].copy()
    errors_cannot_fix.drop(columns='combined_column', inplace=True)

errors_cannot_fix.to_csv('results/allele_cannot_fix_format1.tsv', sep='\t', index=False)


# Separate into error types
error_cols = {'pattern_error', 'invalid_error', 'sequence_error'}
id_cols = set(errors_cannot_fix.columns) - error_cols

errors_cannot_fix = errors_cannot_fix.melt(id_vars=id_cols, value_vars=error_cols, var_name='error_type', value_name='error_info')
errors_cannot_fix = errors_cannot_fix[errors_cannot_fix['error_info'] != '']
errors_cannot_fix = errors_cannot_fix[['systematic_id', 'allele_name', 'allele_description', 'error_type', 'error_info']].copy()

sequence_errors = errors_cannot_fix[errors_cannot_fix.error_type == 'sequence_error'].drop(columns=['error_type']).rename(columns={'error_info': 'sequence_error'})
sequence_errors.to_csv('results/allele_cannot_fix_sequence_errors.tsv', sep='\t', index=False)

errors_cannot_fix[errors_cannot_fix.error_type != 'sequence_error'].to_csv('results/allele_cannot_fix_other_errors.tsv', sep='\t', index=False)
