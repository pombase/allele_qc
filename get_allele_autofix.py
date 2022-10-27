"""
From the output file of 'fix_coordinates' or 'perform_qc' generate two output files:

- --output_auto_fix: the alleles that can be safely renamed (without extra supervision) with the following columns:

    fix_type	systematic_id	allele_name	allele_type	allele_description	change_type_to	rename_to

    - 'fix_type' contains the nature of the error
    - for each allele, 'rename_to' contains the value of 'rename_to' in the input file, unless the coordinates
    have been changed, in which case it contains the value of 'after_coords_rename_to'.

- --output_needs_supervision: the alleles that have sequence errors so may need human double-checking (this file is
    only created if the --input file has the after_coords columns)

"""
import pandas
import argparse
from refinement_functions import seq_error_rename_to

parser = argparse.ArgumentParser(description='Format alleles for auto-fix')
parser.add_argument('--input', default='results/allele_results_after_coordinates.tsv')
parser.add_argument('--output_auto_fix', default='results/allele_auto_fix.tsv')
parser.add_argument('--output_needs_supervision', default='results/allele_needs_supervision.tsv')
args = parser.parse_args()

# Only the alleles that needed fixing
data = pandas.read_csv(args.input, delimiter='\t', na_filter=False)

# Special case for unknown alleles with description unknown:
unknown_in_description = data['allele_description'] == 'unknown'
data.loc[unknown_in_description, 'needs_fixing'] = True

output_data = data[data['needs_fixing'] == True].copy()

# Create a new column to indicate the type of fix
output_data['fix_type'] = ''

# Set the value of fix_type for different types of errors
seq_error = output_data['sequence_error'] != ''
syntax_error = output_data['rename_to'] != ''
type_error = output_data['change_type_to'] != ''

output_data.loc[~seq_error & syntax_error, 'fix_type'] = 'syntax_error'
output_data.loc[~seq_error & (output_data['allele_description'] == 'unknown'), 'fix_type'] = 'unknown_in_description'
output_data.loc[~seq_error & type_error, 'fix_type'] = 'type_error'
output_data.loc[~seq_error & type_error & syntax_error, 'fix_type'] = 'type_error_and_syntax_error'

# Same thing, including columns related to sequence coordinate fixing
if 'after_coords_rename_to' in output_data:
    fixed_after_coord_change = seq_error & (output_data['after_coords_sequence_error'] == '') & (output_data['after_coords_rename_to'] != '')

    output_data.loc[fixed_after_coord_change, 'fix_type'] = 'update_coordinates'
    output_data.loc[fixed_after_coord_change, 'rename_to'] = output_data.loc[fixed_after_coord_change, 'after_coords_rename_to']

    # Partial aminoacid deletions that did not cause errors before or after the change in coordinates
    coordinate_change_no_mutation = ~seq_error & (output_data['after_coords_sequence_error'] == '') & (output_data['after_coords_rename_to'] != '')
    output_data.loc[coordinate_change_no_mutation, 'fix_type'] = 'coordinate_change_no_mutation'

output_data = output_data.sort_values(['fix_type', 'systematic_id', 'allele_name'])

auto_fix = output_data[(output_data['fix_type'] != '') & (output_data['fix_type'] != 'coordinate_change_no_mutation')]
auto_fix[['fix_type', 'systematic_id', 'allele_name', 'allele_type', 'allele_description', 'change_type_to', 'rename_to']].to_csv(args.output_auto_fix, sep='\t', index=False)


if 'after_coords_rename_to' in output_data:

    # When coordinates are off by one
    shift_coordinates_fix = seq_error & ~fixed_after_coord_change & (output_data['after_coords_sequence_error'].str.contains('>') | output_data['sequence_error'].str.contains('>')) & ~output_data['allele_type'].str.contains('nucl')
    shift_coordinates_data = output_data[shift_coordinates_fix].copy()
    shift_coordinates_data['proposed_fix'] = ''
    shift_coordinates_data['after_coords_proposed_fix'] = ''
    for i, row in shift_coordinates_data.iterrows():
        name = row['rename_to'] if row['rename_to'] else row['allele_description']
        shift_coordinates_data.at[i, 'proposed_fix'] = seq_error_rename_to(name, row['sequence_error'])
        # When coordinates are odd by one and updated
        if row['after_coords_sequence_error']:
            shift_coordinates_data.at[i, 'after_coords_proposed_fix'] = seq_error_rename_to(name, row['after_coords_sequence_error'])
        pass
    shift_coordinates_data['fix_type'] = 'shift_coordinates'
    shift_coordinates_data.loc[output_data['after_coords_sequence_error'] != '', 'fix_type'] = 'update_and_shift_coordinates'

    # Append the coordinate_change_no_mutation ones
    coordinate_no_mutation_data = output_data[coordinate_change_no_mutation].copy()
    coordinate_no_mutation_data['proposed_fix'] = coordinate_no_mutation_data['rename_to']
    coordinate_no_mutation_data['after_coords_proposed_fix'] = coordinate_no_mutation_data['after_coords_rename_to']

    needs_supervision_data = pandas.concat([coordinate_no_mutation_data, shift_coordinates_data])
    needs_supervision_data[['fix_type', 'systematic_id', 'allele_name', 'allele_type', 'allele_description', 'change_type_to', 'rename_to', 'after_coords_rename_to', 'sequence_error', 'after_coords_sequence_error', 'proposed_fix', 'after_coords_proposed_fix']].to_csv(args.output_needs_supervision, sep='\t', index=False)
