import pandas
from allele_fixes import multi_shift_fix, old_coords_fix
from allele_fixes import extract_groups_from_targets, shift_coordinates_by_x


def apply_multi_shift_fix(row, genome, target_column):
    # We use at least 4 for the multi-shift
    if (row[target_column].count(',') < 3):
        return ''
    if row['systematic_id'] not in genome:
        return ''
    peptide_seq = genome[row['systematic_id']]['peptide']
    return '|'.join(multi_shift_fix(peptide_seq, row[target_column].split(',')))


def apply_old_coords_fix(row, coordinate_changes_dict, target_column):

    if row['systematic_id'] not in coordinate_changes_dict:
        return '', '', ''
    result = old_coords_fix(coordinate_changes_dict[row['systematic_id']], row[target_column].split(','))
    valid_solutions = pandas.DataFrame([r for r in result if all(r['matches'])])
    if valid_solutions.empty:
        return '', '', ''
    valid_solutions.loc[:, 'values'] = valid_solutions['values'].apply(','.join)
    # There could be several solutions
    return tuple('|'.join(valid_solutions[key]) for key in ['values', 'revision', 'location'])


# These are histone proteins that typically did not count the methionine
histones = ['SPBC1105.11c', 'SPBC1105.12', 'SPAC1834.03c', 'SPAC1834.04', 'SPAC19G12.06c', 'SPBC8D2.03c', 'SPBC8D2.04', 'SPCC622.08c', 'SPCC622.09']


def apply_histone_fix(row, genome, target_column):
    if not (row['systematic_id'] in histones):
        return ''

    groups = extract_groups_from_targets(row[target_column].split(','))
    peptide_seq = genome[row['systematic_id']]['peptide']

    new_coords, shifted_coords_match = shift_coordinates_by_x(groups, peptide_seq, 'number', 1)
    if all(shifted_coords_match):
        return ','.join(new_coords)

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
