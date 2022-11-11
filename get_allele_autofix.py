"""
From the output file of 'fix_coordinates' or 'perform_qc' generate two output files:

- --output_auto_fix: the alleles that can be safely renamed (without extra supervision) with the following columns:

    fix_type	systematic_id	allele_name	allele_type	allele_description	change_type_to	change_description_to

    - 'fix_type' contains the nature of the error
    - for each allele, 'change_description_to' contains the value of 'change_description_to' in the input file, unless the coordinates
    have been changed, in which case it contains the value of 'after_coords_change_description_to'.

- --output_needs_supervision: the alleles that have sequence errors so may need human double-checking (this file is
    only created if the --input file has the after_coords columns)

- --output_cannot_fix: the alleles that cannot be fixed with this pipeline.

"""
import pandas
import argparse
from refinement_functions import seq_error_change_description_to

parser = argparse.ArgumentParser(description='Format alleles for auto-fix')
parser.add_argument('--input', default='results/allele_results_after_coordinates.tsv')
parser.add_argument('--output_auto_fix', default='results/allele_auto_fix.tsv')
parser.add_argument('--output_needs_supervision', default='results/allele_needs_supervision.tsv')
parser.add_argument('--output_cannot_fix', default='results/allele_cannot_fix.tsv')
args = parser.parse_args()

# Only the alleles that needed fixing
data = pandas.read_csv(args.input, delimiter='\t', na_filter=False)

# Special case for unknown alleles with description unknown:
unknown_in_description = data['allele_description'] == 'unknown'
data.loc[unknown_in_description, 'needs_fixing'] = True

allele_error_data = data[data['needs_fixing'] == True].copy()

# Create a new column to indicate the type of fix
allele_error_data['fix_type'] = ''

# Extra column that is used in manual fixing sometimes
allele_error_data['change_name_to'] = ''

# Set the value of fix_type for different types of errors
seq_error = allele_error_data['sequence_error'] != ''
syntax_error = allele_error_data['change_description_to'] != ''
type_error = allele_error_data['change_type_to'] != ''

allele_error_data.loc[~seq_error & syntax_error, 'fix_type'] = 'syntax_error'
allele_error_data.loc[~seq_error & (allele_error_data['allele_description'] == 'unknown'), 'fix_type'] = 'unknown_in_description'
allele_error_data.loc[~seq_error & type_error, 'fix_type'] = 'type_error'
allele_error_data.loc[~seq_error & type_error & syntax_error, 'fix_type'] = 'type_error_and_syntax_error'

# Same thing, including columns related to sequence coordinate fixing
if 'after_coords_change_description_to' in allele_error_data:
    fixed_after_coord_change = seq_error & (allele_error_data['after_coords_sequence_error'] == '') & (allele_error_data['after_coords_change_description_to'] != '')

    allele_error_data.loc[fixed_after_coord_change, 'fix_type'] = 'update_coordinates'
    allele_error_data.loc[fixed_after_coord_change, 'change_description_to'] = allele_error_data.loc[fixed_after_coord_change, 'after_coords_change_description_to']

    # Partial aminoacid deletions that did not cause errors before or after the change in coordinates
    coordinate_change_no_mutation = ~seq_error & (allele_error_data['after_coords_sequence_error'] == '') & (allele_error_data['after_coords_change_description_to'] != '')
    allele_error_data.loc[coordinate_change_no_mutation, 'fix_type'] = 'coordinate_change_no_mutation'

    # Changing the reference sequence causes sequence errors and before there were no errors
    coordinate_change_creates_error = ~seq_error & (allele_error_data['after_coords_sequence_error'] != '')
    allele_error_data.loc[coordinate_change_creates_error, 'fix_type'] = 'coordinate_change_creates_error'

allele_error_data = allele_error_data.sort_values(['fix_type', 'systematic_id', 'allele_name'])
auto_fix = (allele_error_data['fix_type'] != '') & ~coordinate_change_no_mutation & ~coordinate_change_creates_error
auto_fix_data = allele_error_data.loc[auto_fix]
auto_fix_data[['fix_type', 'systematic_id', 'allele_name', 'reference', 'allele_type', 'allele_description', 'change_type_to', 'change_description_to', 'change_name_to']].to_csv(args.output_auto_fix, sep='\t', index=False)

if 'after_coords_change_description_to' not in allele_error_data:
    allele_error_data[~auto_fix].to_csv(args.output_cannot_fix, sep='\t', index=False)
    exit()

# When coordinates are off by one
shift_coordinates_fix = seq_error & ~fixed_after_coord_change & (allele_error_data['after_coords_sequence_error'].str.contains('>') | allele_error_data['sequence_error'].str.contains('>')) & ~allele_error_data['allele_type'].str.contains('nucl')
shift_coordinates_data = allele_error_data.loc[shift_coordinates_fix].copy()
shift_coordinates_data['proposed_fix'] = ''
shift_coordinates_data['after_coords_proposed_fix'] = ''
for i, row in shift_coordinates_data.iterrows():
    name = row['change_description_to'] if row['change_description_to'] else row['allele_description']
    shift_coordinates_data.at[i, 'proposed_fix'] = seq_error_change_description_to(name, row['sequence_error'])
    # When coordinates are odd by one and updated
    if row['after_coords_sequence_error']:
        shift_coordinates_data.at[i, 'after_coords_proposed_fix'] = seq_error_change_description_to(name, row['after_coords_sequence_error'])
    pass
shift_coordinates_data['fix_type'] = 'shift_coordinates'
shift_coordinates_data.loc[allele_error_data['after_coords_sequence_error'] != '', 'fix_type'] = 'update_and_shift_coordinates'

# When changing the reference sequence causes sequence errors and before there were no errors
coordinate_change_creates_error_data = allele_error_data.loc[coordinate_change_creates_error].copy()
coordinate_change_creates_error_data['proposed_fix'] = ''
coordinate_change_creates_error_data['after_coords_proposed_fix'] = ''


# Append the coordinate_change_no_mutation ones
coordinate_no_mutation_data = allele_error_data.loc[coordinate_change_no_mutation].copy()
coordinate_no_mutation_data['proposed_fix'] = ''
coordinate_no_mutation_data['after_coords_proposed_fix'] = ''

needs_supervision_data = pandas.concat([coordinate_no_mutation_data, coordinate_change_creates_error_data, shift_coordinates_data])
needs_supervision_data[['fix_type', 'systematic_id', 'allele_name', 'reference', 'allele_type', 'allele_description', 'change_type_to', 'change_description_to', 'change_name_to', 'after_coords_change_description_to', 'sequence_error', 'after_coords_sequence_error', 'proposed_fix', 'after_coords_proposed_fix']].to_csv(args.output_needs_supervision, sep='\t', index=False)

allele_error_data.loc[~auto_fix & ~shift_coordinates_fix & ~coordinate_change_no_mutation & ~coordinate_change_creates_error].to_csv(args.output_cannot_fix, sep='\t', index=False)
