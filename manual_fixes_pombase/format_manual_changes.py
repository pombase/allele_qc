import pandas

autofix_data = pandas.read_csv('../results/allele_auto_fix.tsv', sep='\t', na_filter=False)
manual_data1 = pandas.read_excel('allele_cannot_fix.xlsx', sheet_name=None, na_filter=False)
manual_data2 = pandas.read_excel('allele_needs_supervision.xlsx', sheet_name=None, na_filter=False)

manual_data = pandas.concat(list(manual_data1.values()) + list(manual_data2.values()))

cols2keep = ['systematic_id', 'allele_name', 'manual_fix_allele_description', 'manual_fix_allele_name', 'manual_fix_allele_type', 'comment']

manual_data = manual_data[cols2keep].copy()
manual_data.rename(columns={'comment': 'manual_fix_comment'})

# Alleles that have to be removed (added to separate file)
alleles2remove = manual_data.manual_fix_allele_description.str.contains('remove')
manual_data[alleles2remove].to_csv('manual_alleles2remove.tsv', sep='\t', index=False)
manual_data = manual_data[~alleles2remove].copy()

# Those without fix
no_fix = manual_data.manual_fix_allele_description.str.contains('cannot fix') | \
    manual_data[['manual_fix_allele_description', 'manual_fix_allele_name', 'manual_fix_allele_type', 'comment']].eq('').all(axis=1)

manual_data = manual_data[~no_fix].copy()

# Load the autofix, and look for conflicts
merged_data = autofix_data.merge(manual_data[manual_data.manual_fix_allele_description != ''], on='allele_name', how='inner')
merged_data.fillna('', inplace=True)
conflict_logi = merged_data.manual_fix_allele_description != merged_data.change_description_to
conflict_data = merged_data.loc[conflict_logi, ['systematic_id_x', 'allele_description', 'manual_fix_allele_description', 'solution_index', 'change_description_to', 'auto_fix_comment']]
if not conflict_data.empty:
    print('\033[1;33mConflict between manual and auto fix data printed to manual_auto_fix_conflict_data.tsv\033[0m')
    conflict_data.to_csv('manual_auto_fix_conflict_data.tsv', sep='\t', index=False)


# In some cases I had not added the systematic_id, so we fill from the original allele data
original_allele_data = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)
manual_data.fillna('', inplace=True)
missing_systematic_id = manual_data['systematic_id'] == ''
manual_data.loc[missing_systematic_id, 'systematic_id'] = manual_data.loc[missing_systematic_id, 'allele_name'].apply(lambda name: original_allele_data.loc[original_allele_data.allele_name == name].iloc[0])

manual_data.to_csv('manual_changes_formatted.tsv', sep='\t', index=False)
