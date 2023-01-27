"""
Updates the modification files in the svn repository with the current
"""

import os
import glob
import pandas
import tempfile

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
all_applied_fixes = list()

# Temporary directory to store the files that will eventually replace the current ones
temp_dir = tempfile.gettempdir()

for modification_file in glob.glob(modification_folder + '/PMID*.tsv'):
    file_name = os.path.split(modification_file)[-1]
    print(file_name)
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
    fixes_for_this_paper = fixes.loc[(fixes.reference == pmid) & (fixes.date == data['date'][0]), :].copy()

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

    # We store the file in a temp directory (we want to replace the original
    # files only if all checks have passed).
    data_fixed.to_csv(os.path.join(temp_dir, file_name), sep='\t', index=False)


for f in glob.glob(os.path.join(temp_dir, '*.tsv')):
    file_name = os.path.basename(f)
    os.rename(f, os.path.join(modification_folder, file_name))
    print(file_name, 'overwritten')
