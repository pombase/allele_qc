import pandas
from allele_fixes import multi_shift_fix, old_coords_fix, shift_coordinates_by_x, position_or_index_exists
import re


def apply_multi_shift_fix(row, genome, target_column):
    # We use at least 4 for the multi-shift, and they must contain the expected
    # residue (e.g. A123, A123V, but not 12-14).
    all_targets_with_residue = [i for i in row[target_column].split(',') if i[0].isalpha()]
    # Also they should be referring to unique sites (e.g. not A123V, A123P, A123L, ect.)
    unique_residue_references = set(i[:-1] if i[-1].isalpha() else i for i in all_targets_with_residue)
    if len(unique_residue_references) < 4:
        return ''
    if row['systematic_id'] not in genome:
        return ''
    if 'peptide' not in genome[row['systematic_id']]:
        print(f'no peptide found for {row["systematic_id"]}')
        return ''

    peptide_seq = genome[row['systematic_id']]['peptide']
    return '|'.join(multi_shift_fix(peptide_seq, row[target_column].split(',')))


def apply_old_coords_fix(row, coordinate_changes_dict, target_column):

    if row['systematic_id'] not in coordinate_changes_dict:
        return '', '', ''
    possible_fixes = old_coords_fix(coordinate_changes_dict[row['systematic_id']], row[target_column].split(','))

    if len(possible_fixes) == 0:
        return '', '', ''
    # We convert a list of dictionaries to dataframe to concatenate all values, all revisions, all location using |
    possible_fixes = pandas.DataFrame(possible_fixes)
    # There could be several solutions
    return tuple('|'.join(possible_fixes[key]) for key in ['values', 'revision', 'location'])


# These are histone proteins that typically did not count the methionine
histones = ['SPBC1105.11c', 'SPBC1105.12', 'SPAC1834.03c', 'SPAC1834.04', 'SPAC19G12.06c', 'SPBC8D2.03c', 'SPBC8D2.04', 'SPCC622.08c', 'SPCC622.09', 'SPBC11B10.10c', 'SPBC1105.17']


def apply_histone_fix(row, genome, target_column):
    if not (row['systematic_id'] in histones):
        return ''

    shifted_targets = shift_coordinates_by_x(row[target_column], 1)
    peptide_seq = genome[row['systematic_id']]['peptide']

    if all(position_or_index_exists(shifted_target, peptide_seq) for shifted_target in shifted_targets.split(',')):
        return shifted_targets

    return ''


def get_preferred_fix(row):
    """
    Based on a hierarchy
    """

    if row['histone_fix']:
        return row['histone_fix'], 'histone_fix'
    if row['old_coords_fix']:
        fix_reasons = list()
        for revision, location in zip(row["old_coords_revision"].split('|'), row["old_coords_location"].split('|')):
            fix_reasons.append(f'old_coords_fix, revision {revision}: {location}')
        return row['old_coords_fix'], '|'.join(fix_reasons)
    if row['multi_shift_fix']:
        return row['multi_shift_fix'], 'multi_shift_fix'

    return '', ''


def format_auto_fix(row, target_column, syntax_error_column):

    syntax_error = row[syntax_error_column] != ''
    sequence_position = row[target_column] if not syntax_error else row[syntax_error_column]

    # Cases where there was no error in the sequence, but there might be a syntax error that needs fixing
    if not row['auto_fix_to']:
        return (row[syntax_error_column], 'syntax_error') if syntax_error else ('', '')

    # There are cases where multiple fixes could be possible, and cases where changes were made
    # and reverted several times, so the same coordinate (or updated depending on genome sequence)
    # may exist
    possible_fixes = list()
    for all_change_to in row['auto_fix_to'].split('|'):
        this_fix = sequence_position
        for change_from, change_to in zip(row['auto_fix_from'].split(','), all_change_to.split(',')):
            this_fix = this_fix.replace(change_from, change_to)
        possible_fixes.append(this_fix)

    # Case 1 -> all old coordinates give the same residue
    if len(set(possible_fixes)) == 1:
        return possible_fixes[0], '/'.join(row['auto_fix_comment'].split('|'))

    # Case 2 -> different old coordinates give different residues
    else:
        return '|'.join(possible_fixes), row['auto_fix_comment']


def print_warnings(data: pandas.DataFrame):

    print('\033[0;32mmixed case\033[0m')
    print()
    for i, row in data.iterrows():
        desc = row['allele_description']
        # Combinations of upper and lower case might have been weird patterns
        if re.search('[a-z]', desc) and re.search('[A-Z]', desc):
            new_desc = row['change_description_to']
            if new_desc == '':
                print(desc)
            else:
                print(desc, '    >>>>>    ', new_desc)

    print()
    print('\033[0;32mnt or aa\033[0m')
    print()
    for i, row in data.iterrows():
        desc = row['allele_description']

        # Presence of aa or nt might mean nucleotides or aminoacids
        if 'aa' in desc or 'nt' in desc:
            new_desc = row['change_description_to']
            if new_desc == '':
                print(desc)
            else:
                print(desc, '    >>>>>    ', new_desc)


def apply_name_fix(row):

    if row['sequence_error'] == '' or row['auto_fix_comment'] == 'histone_fix':
        return ''

    new_description = row['change_description_to']
    old_description = row['allele_description']
    name = row['allele_name']
    if old_description in name:
        return name.replace(old_description, new_description)
    return ''
