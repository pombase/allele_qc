"""
Reformat the manual changes to match the allele_auto_fix + extra column with a comment

fix_type	systematic_id	allele_name	reference	allele_type	allele_description	change_type_to	change_description_to	change_name_to  comment

Also create a file with the ones that need to be deleted

"""

import pandas
from numpy import NaN
input_alleles = pandas.read_csv('../data/alleles.tsv', sep='\t', na_filter=False)

identified_alleles = pandas.read_excel('allele_cannot_fix.xlsx', sheet_name='identified', na_filter=False)

identified_alleles['fix_type'] = 'manual_fix'

# Drop old fixes that still resulted in errors
identified_alleles = identified_alleles.drop(columns=['change_description_to', 'change_type_to'])

unidentified_alleles = pandas.read_excel('allele_cannot_fix.xlsx', sheet_name='not identified')

# Fill columns from the original alleles
for col_name in ['systematic_id', 'reference', 'allele_type', 'allele_description']:
    unidentified_alleles[col_name] = NaN
unidentified_alleles = unidentified_alleles.fillna(input_alleles)

unidentified_alleles['fix_type'] = 'unnoticed'

manual_data = pandas.concat([identified_alleles, unidentified_alleles])

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
changes_applied['allele_name'][changes_applied['change_name_to'] != ''] = changes_applied['change_name_to']
changes_applied['allele_description'][changes_applied['change_description_to'] != ''] = changes_applied['change_description_to']

changes_applied = changes_applied.drop(columns=['change_type_to', 'change_description_to', 'change_name_to'])
changes_applied.to_csv('manual_changes_applied.tsv', sep='\t', index=False)
