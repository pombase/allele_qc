import pandas

# Original list of alleles, used to verify that no new allele names have been included
input_alleles = pandas.read_csv('data/alleles.tsv', sep='\t', na_filter=False)

# All the list of manual data (not all of them have been corrected)
manual_data = pandas.read_csv('manual_fixes_pombase/manual_changes_formatted.tsv', sep='\t', na_filter=False)

# Verify that no new names have been added
new_names = ~manual_data['allele_name'].isin(input_alleles['allele_name'])
if any(new_names):
    raise ValueError('New alleles exist', list(manual_data['allele_name'][new_names]))

# Quality control of the manual data
manual_data_qc = pandas.read_csv('manual_fixes_pombase/allele_results_manual_changes.tsv', sep='\t', na_filter=False)

# Alleles that can be auto-fixed
autofix_data = pandas.read_csv('results/allele_auto_fix.tsv', sep='\t', na_filter=False)

# Verify that no new names have been added
new_names = ~autofix_data['allele_name'].isin(input_alleles['allele_name'])
if any(new_names):
    raise ValueError('New alleles exist', list(autofix_data['allele_name'][new_names]))

# Format columns
autofix_data['comment'] = ''
autofix_data = autofix_data[['fix_type', 'systematic_id', 'allele_name', 'reference', 'allele_type', 'allele_description', 'change_type_to', 'change_description_to', 'change_name_to', 'comment']].copy()

# The manual fixes that passe the qc and therefore can be used in batch-fix
manual_data_correct = manual_data[manual_data.allele_name.isin(manual_data_qc.allele_name)].copy()

# remove the rows that have a manual fix from the autofix (some are overwritten)
manually_fixed = autofix_data['allele_name'].isin(manual_data_correct['allele_name'])
autofix_data = autofix_data.loc[~manually_fixed].copy()

# Concatenate both to create the batch-fix file
batch_fix = pandas.concat([autofix_data, manual_data_correct])
batch_fix.to_csv('output_files/batch_fix.tsv', sep='\t', index=False)

# Report a file with unfixed sequence errors

# Load the file with all errors

error_data = pandas.read_csv('results/allele_results_errors.tsv', sep='\t', na_filter=False)

# See which ones have not been fixed
unfixed_sequence_errors = (error_data.sequence_error != '') & ~error_data.allele_name.isin(batch_fix.allele_name)
error_data[unfixed_sequence_errors].to_csv('output_files/unfixed_sequence_errors.tsv', sep='\t', index=False)
