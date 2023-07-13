import pandas

existing_fixes = pandas.read_csv('change_log/protein_modification_auto_fix_external_data_21032023.tsv', sep='\t', na_filter=False)
new_fixes = pandas.read_csv('results/protein_modification_auto_fix.tsv', sep='\t', na_filter=False)

# drop the row where systematic_id is SPAC1486.05
row2remove = (new_fixes['systematic_id'] == 'SPAC1486.05') & (new_fixes['change_sequence_position_to'] == 'K1001')
new_fixes = new_fixes[~row2remove].copy()

new_fixes.drop(columns=['solution_index'], inplace=True)

# We drop sequence error columns in both, because the content changed
# between versions (we use no ">" to indicate close residues that match anymore)
new_fixes.drop(columns=['sequence_error'], inplace=True)
existing_fixes.drop(columns=['sequence_error'], inplace=True)


columns4merge = list(new_fixes.columns)
columns4merge.remove('change_sequence_position_to')
result = existing_fixes.merge(new_fixes, on=columns4merge, how='left', indicator=True)

# First of all, we make sure that the only column in which they differ is change_sequence_position_to
print('are they different in any other column?', not (result[result['_merge'] == 'left_only'].empty))

# Now drop the columns where the change_position_to is the same
result = result[result['change_sequence_position_to_x'] != result['change_sequence_position_to_y']].copy()

result.to_csv('different.tsv', sep='\t', index=False)
