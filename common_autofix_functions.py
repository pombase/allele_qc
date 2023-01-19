import pandas
from allele_fixes import multi_shift_fix, old_coords_fix, shift_coordinates_by_x, position_or_index_exists


def apply_multi_shift_fix(row, genome, target_column):
    # We use at least 4 for the multi-shift
    if (row[target_column].count(',') < 3):
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
histones = ['SPBC1105.11c', 'SPBC1105.12', 'SPAC1834.03c', 'SPAC1834.04', 'SPAC19G12.06c', 'SPBC8D2.03c', 'SPBC8D2.04', 'SPCC622.08c', 'SPCC622.09', 'SPBC11B10.10c']


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

    # There are edge cases where multiple fixes could be possible
    possible_fixes = set()
    for all_change_to in row['auto_fix_to'].split('|'):
        this_fix = sequence_position
        for change_from, change_to in zip(row['auto_fix_from'].split(','), all_change_to.split(',')):
            this_fix = this_fix.replace(change_from, change_to)
        possible_fixes.add(this_fix)

    return '|'.join(possible_fixes), row['auto_fix_comment']
