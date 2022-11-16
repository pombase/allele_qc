"""
Reformat the manual changes to match the allele_auto_fix + extra column with a comment

fix_type	systematic_id	allele_name	reference	allele_type	allele_description	change_type_to	change_description_to	change_name_to  comment

Also create a file with the ones that need to be deleted

"""

import pandas

input_alleles = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)

identified_alleles = pandas.read_excel('allele_cannot_fix.xlsx', sheet_name='identified', na_filter=False)

identified_alleles['fix_type'] = 'manual_fix_other'
identified_alleles.loc[identified_alleles['invalid_error'] != '', 'fix_type'] = 'manual_fix_invalid_error'
identified_alleles.loc[identified_alleles['sequence_error'] != '', 'fix_type'] = 'manual_fix_sequence_error'

unidentified_alleles = pandas.read_excel('allele_cannot_fix.xlsx', sheet_name='not identified')

# Fill columns from the original alleles
cols = ['systematic_id', 'reference', 'allele_type', 'allele_description']
for col in cols:
    unidentified_alleles[col] = ''

for i, row in unidentified_alleles.iterrows():
    for col in cols:
        unidentified_alleles.at[i, col] = list(input_alleles[input_alleles.allele_name == row.allele_name][col])[0]

unidentified_alleles['fix_type'] = 'unnoticed'

# Formatting here is required
need_supervision = pandas.read_excel('allele_needs_supervision.xlsx', na_filter=False)

# Remove the ones that were either fixed in canto or can't be fixed
need_supervision = need_supervision[~need_supervision.manual_fix_allele_description.isin(['skip', 'asked'])]

manual_data = pandas.concat([identified_alleles, unidentified_alleles, need_supervision])

# Drop old fixes that still resulted in errors
manual_data = manual_data.drop(columns=['change_description_to', 'change_type_to'])

manual_data.rename(columns={
    'manual_fix_allele_description': 'change_description_to',
    'manual_fix_allele_name': 'change_name_to',
    'manual_fix_allele_type': 'change_type_to',
}, inplace=True)

manual_data = manual_data.fillna('')

manual_data = manual_data[['fix_type', 'systematic_id', 'allele_name', 'reference', 'allele_type', 'allele_description', 'change_type_to', 'change_description_to', 'change_name_to', 'comment']]

manual_data.to_csv('manual_changes_formatted.tsv', sep='\t', index=False)

changes_applied = manual_data.copy()

changes_applied['allele_type'][changes_applied['change_type_to'] != ''] = changes_applied['change_type_to']
changes_applied['allele_description'][changes_applied['change_description_to'] != ''] = changes_applied['change_description_to']

# I don't do this, because we want to preserve the original names to see which one gave errors.
# changes_applied['allele_name'][changes_applied['change_name_to'] != ''] = changes_applied['change_name_to']

changes_applied = changes_applied.drop(columns=['change_type_to', 'change_description_to', 'change_name_to'])
changes_applied.to_csv('manual_changes_applied.tsv', sep='\t', index=False)
