"""
Usage:

python update_svn_modification_files.py svn_repo_path

svn_repo_path is the root directory of pombe svn (where branches are trunk dirs are).

The script applies the changes in results/protein_modification_auto_fix.tsv to
the modification files in the svn repository (trunk/external_data/modification_files).
The files with errors corrected are saved to a folder in this repo called
svn_protein_modification_files_corrected (included in .gigtignore).

In the same folder, along with the corrected files, a file called fixes_applied.tsv
is created with the lines in results/protein_modification_auto_fix.tsv
that were used.

Then you can simply replace the svn files by running:
mv svn_protein_modification_files_corrected/PMID*.tsv svn_repo_path/trunk/external_data/modification_files/
"""

import os
import glob
import pandas
import numpy as np
import warnings

# Column names are not consistent across files
column_names = [
    'systematic_id',
    'primary_name',
    'modification',
    'evidence',
    'sequence_position',
    'annotation_extension',
    'reference',
    'taxon',
    'date'
]

svn_folder = '/Users/manu/Documents/OpenSource/pombe-embl'
modification_sub_path = 'trunk/external_data/modification_files'

modification_folder = os.path.join(svn_folder, modification_sub_path)

fixes = pandas.read_csv('results/protein_modification_auto_fix.tsv', sep='\t', na_filter=False)
fixed_references = set(fixes.reference)

# This is an extra control, to make sure we are not replacing by empty strings (no solution found),
# Empty strings should never be included in protein_modification_auto_fix.tsv,
# and be in protein_modification_cannot_fix.tsv, but just in case
if any(fixes.change_sequence_position_to == ''):
    raise ValueError(f'Error in {fixes}: empty string in change_sequence_position_to')

# Create output folder if it does not exist:
output_folder = 'svn_protein_modification_files_corrected'
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

# Logical indexes of the rows in fixes that are used, updated at each iteration
all_applied_fixes = np.zeros_like(fixes['systematic_id'], dtype=bool)

# Logical indexes of rows that contain several solutions
multi_solution_fixes = fixes['solution_index'] != ''
fixes.drop(columns='solution_index', inplace=True)
multi_solution_skipped_warnings = list()

for modification_file in glob.glob(modification_folder + '/PMID*.tsv'):
    file_name = os.path.split(modification_file)[-1]
    print(file_name)

    # Read the comments of the file to write them again later
    comments = list()
    with open(modification_file) as ins:
        for line in ins:
            if line.startswith('#'):
                comments.append(line)
            else:
                break

    data = pandas.read_csv(modification_file, sep='\t', na_filter=False, comment='#')

    original_column_names = list(data.columns)
    data.columns = column_names

    # For some reason the merge below changes the order of rows, so we keep the value here
    data['original_index'] = data.index

    # There should only be one PMID per file
    if (len(pandas.unique(data['reference'])) != 1) or (len(pandas.unique(data['date'])) != 1):
        raise ValueError(f'Error in {file_name}: more than one PMID or date')

    pmid = pandas.unique(data['reference'])[0]

    if pmid not in fixed_references:
        continue

    print('  errors found')

    # Important to include also the date, because some of the changes were made in canto, so they belong
    # to the same PMID, but not in this file.
    rows_for_this_paper = (fixes.reference == pmid) & (fixes.date == data['date'][0])
    # Check if there are multi-solution for this paper, and skip those rows:
    if any(multi_solution_fixes & rows_for_this_paper):
        multi_solution_skipped_warnings.append(file_name)

    rows_for_this_paper = rows_for_this_paper & ~multi_solution_fixes
    fixes_for_this_paper = fixes.loc[rows_for_this_paper, :].copy()
    all_applied_fixes = all_applied_fixes | rows_for_this_paper

    # We create a column combining the systematic_id + residue, which uniquely identifies the position
    # in the genome.
    fixes_for_this_paper['combined_column'] = fixes_for_this_paper.apply(lambda r: r['systematic_id'] + '|' + r['sequence_position'], axis=1)
    data['combined_column'] = data.apply(lambda r: r['systematic_id'] + '|' + r['sequence_position'], axis=1)

    # We verify that none of those combinations in auto_fix is missing from the dataset.
    in_fixes_only = ~fixes_for_this_paper['combined_column'].isin(set(data['combined_column']))
    if any(in_fixes_only):
        raise ValueError(f'Error: the following combinations of residue + systematic id are only present in protein_modification_auto_fix.tsv\n{fixes_for_this_paper.loc[in_fixes_only,"combined_column"]}')
    fixes_for_this_paper.drop(columns='combined_column', inplace=True)
    data.drop(columns='combined_column', inplace=True)

    # We replace by the right values
    data_fixed = data.merge(
        fixes_for_this_paper[['systematic_id', 'sequence_position', 'change_sequence_position_to']].drop_duplicates(),
        on=['systematic_id', 'sequence_position'], how='outer'
    )

    # Sort so that the original index of rows is kept and drop the column
    data_fixed.sort_values('original_index', inplace=True)
    data_fixed.drop(columns='original_index', inplace=True)

    # We also drop original_index in data, because later we will check that the merge does not change
    # the shape of the dataset
    data.drop(columns='original_index', inplace=True)

    data_fixed.fillna('', inplace=True)

    data_fixed.loc[:, 'sequence_position'] = data_fixed.loc[:, 'change_sequence_position_to']
    data_fixed.drop(columns='change_sequence_position_to', inplace=True)

    # The resulting file should have the same number of rows as the original one
    if (data_fixed.shape != data.shape):
        raise AssertionError('The fixes change the length of the dataset')

    # Revert to original column names
    data_fixed.columns = original_column_names

    # Write the file with the original comments and column names
    with open(os.path.join(output_folder, file_name), 'w') as ins:
        ins.writelines(comments)
        data_fixed.to_csv(ins, sep='\t', index=False)

# Write to a file the changes in fixes that have been applied
fixes[all_applied_fixes].to_csv(os.path.join(output_folder, 'fixes_applied.tsv'), sep='\t', index=False)

# Print a warning if any row with multiple solutions was skipped
for f in multi_solution_skipped_warnings:
    warnings.warn(f'\033[1;33mWarning: Some lines in {f} skipped, contained multiple solutions. You have to resolve the conflict by hand \033[0m')
