"""
Reformat the manual changes to match the allele_auto_fix + extra column with a comment

fix_type	systematic_id	allele_name	reference	allele_type	allele_description	change_type_to	change_description_to	change_name_to  comment

Also create a file with the ones that need to be deleted

"""

import pandas

identified_alleles = pandas.read_excel('allele_cannot_fix.xlsx', sheet_name='identified')

identified_alleles['fix_type'] = 'manual_fix'

identified_alleles = identified_alleles.drop(columns=['change_description_to', 'change_type_to'])

identified_alleles.rename(columns={
    'manual_fix_allele_description': 'change_description_to',
    'manual_fix_allele_name': 'change_name_to',
    'manual_fix_allele_type': 'change_type_to',
}, inplace=True)


unidentified_alleles = pandas.read_excel('allele_cannot_fix.xlsx', sheet_name='not identified')

unidentified_alleles.rename(columns={
    'manual_fix_allele_description': 'change_description_to',
}, inplace=True)

unidentified_alleles['fix_type'] = 'unnoticed'

out_data = pandas.concat([identified_alleles, unidentified_alleles])

out_data[['fix_type', 'systematic_id', 'allele_name', 'reference', 'allele_type', 'allele_description', 'change_type_to', 'change_description_to', 'change_name_to', 'comment']].to_csv('../allele_changes_pombase/formatted_manual_changes.tsv', sep='\t', index=False)
