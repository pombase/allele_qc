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

# Create output folder if it does not exist:
output_folder = 'svn_protein_modification_files_corrected'
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

# Logical index of the rows in fixes that are used, updated at each iteration
all_applied_fixes = np.zeros_like(fixes['systematic_id'], dtype=bool)

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
    # We also drop in data, because later we will check that the merge does not change
    # the shape of the dataset
    data.drop(columns='original_index', inplace=True)

    data_fixed.fillna('', inplace=True)
    rows2fix = data_fixed.change_sequence_position_to != ''
    data_fixed.loc[rows2fix, 'sequence_position'] = data_fixed.loc[rows2fix, 'change_sequence_position_to']
    data_fixed.drop(columns='change_sequence_position_to', inplace=True)

    # The resulting file should have the same number of rows as the original one
    if (data_fixed.shape != data.shape):
        raise AssertionError('The fixes change the length of the dataset')

    # Revert to original column names
    data_fixed.columns = original_column_names

    with open(os.path.join(output_folder, file_name), 'w') as ins:
        ins.writelines(comments)
        data_fixed.to_csv(ins, sep='\t', index=False)

fixes[all_applied_fixes].to_csv(os.path.join(output_folder, 'fixes_applied.tsv'), sep='\t', index=False)
