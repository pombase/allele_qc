# %%
import pandas
import argparse

parser = argparse.ArgumentParser(description='Format alleles for auto-fix')
parser.add_argument('--input', default='results/allele_results_after_coordinates.tsv')
parser.add_argument('--output', default='results/allele_auto_fix.tsv')
args = parser.parse_args()

data = pandas.read_csv(args.input, delimiter='\t', na_filter=False)

output_data = data[data['needs_fixing'] == True].copy()

output_data['fix_type'] = ''

seq_error = output_data['sequence_error'] != ''
syntax_error = output_data['rename_to'] != ''
type_error = output_data['change_type_to'] != ''

output_data.loc[~seq_error & syntax_error, 'fix_type'] = 'syntax_error'
output_data.loc[~seq_error & type_error, 'fix_type'] = 'type_error'
output_data.loc[~seq_error & type_error & syntax_error, 'fix_type'] = 'type_error & syntax_error'

fixed_after_coord_change = seq_error & (output_data['after_coords_sequence_error'] == '') & (output_data['after_coords_rename_to'] != '')

output_data.loc[fixed_after_coord_change, 'fix_type'] = 'wrong_coordinates'
output_data.loc[fixed_after_coord_change, 'rename_to'] = output_data.loc[fixed_after_coord_change, 'after_coords_rename_to']

# Drop the ones that can't be fixed
output_data = output_data[output_data['fix_type'] != '']

output_data = output_data.sort_values(['fix_type', 'systematic_id', 'allele_name'])
output_data[['fix_type', 'systematic_id', 'allele_name', 'allele_type', 'allele_description', 'change_type_to', 'rename_to']].to_csv(args.output, sep='\t', index=False)

# Exclude data that still has errors in sequence
