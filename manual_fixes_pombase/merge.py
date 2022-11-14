autofix_data = pandas.read_csv('../results/allele_auto_fix.tsv', sep='\t', na_filter=False)
autofix_data['comment'] = ''
autofix_data = autofix_data[['fix_type', 'systematic_id', 'allele_name', 'reference', 'allele_type', 'allele_description', 'change_type_to', 'change_description_to', 'change_name_to', 'comment']].copy()


# remove the rows that have a manual fix
manually_fixed = autofix_data['allele_description'].isin(manual_data['allele_description'])
autofix_data = autofix_data.loc[~manually_fixed].copy()

output_data = pandas.concat([manual_data, autofix_data])
output_data = output_data.fillna('')

# Verify that no new names have been added
new_names = ~output_data['allele_name'].isin(input_alleles['allele_name'])
if any(new_names):
    raise ValueError('New alleles exist', list(output_data['allele_name'][new_names]))

# Export a file with the changes to be ingested by the database
output_data.to_csv('formatted_changes.tsv', sep='\t', index=False)

# Export a file with the applied changes to verify that the fixes worked (should give no errors)

changes_applied = output_data.copy()

changes_applied['allele_type'][changes_applied['change_type_to'] != ''] = changes_applied['change_type_to']
changes_applied['allele_name'][changes_applied['change_name_to'] != ''] = changes_applied['change_name_to']
changes_applied['allele_description'][changes_applied['change_description_to'] != ''] = changes_applied['change_description_to']

changes_applied.drop(columns=['allele_type', 'allele_name', 'allele_description'])
changes_applied.to_csv('changes_applied.tsv', sep='\t', index=False)
