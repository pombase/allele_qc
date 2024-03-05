import re

aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

ctd_repetition_regex = r'\(r\d+-r\d+(?:-\d+)?\)'
ctd_repetition_regex_with_groups = r'\(r(\d+)-r(\d+)(?:-(\d+))?\)'
ctd_mutation_regex = f'{aa}\d+{aa}(?:{ctd_repetition_regex})?'
ctd_deletion_regex = f'(?:delta|\u0394|∆)(?:{ctd_repetition_regex})?'

ctd_abbreviations = f'(?:CTD_S2|CTD_T4|CTD_S5|CTD_S7)'

# This dictionary has keys that are the systematic IDs of the genes that have CTDs,
# the value is a dictionary with:
# - positions: a dictionary with keys that are the residues of the canonical repeat
#   e.g. YSPTSPS for SPBC28F2.12, and the value is the repeats in which the residue is
#   present, zero-based index.
# - shifts: a tuple of two integers the first value is the position in the protein where
#   an aminoacid is missing (-1) or inserted (+1) with respect to the canonical repeat.
#   There is only one in rpb1, but a list of tuples could be used if more are needed.
# - start: the position of the first residue of the CTD in the protein (zero-based index)
# - length: the length of the CTD

rpb1_dictionary = {
    'positions': {
        'Y1': range(29),
        'S2': [1, 2] + list(range(4, 29)),
        'P3': [1, ] + list(range(3, 8)) + list(range(9, 29)),
        'T4': [0, 2] + list(range(4, 29)),
        'S5': range(29),
        'P6': range(29),
        'S7': [0, 3] + list(range(4, 29)),
    },
    'shift': (1566, -1),
    'start': 1550,
    'length': 202,
    'nb_repeats': 29,
}


def apply_shift(match: re.match, shift: tuple[int, int]):
    num_str = match.group()
    num = int(num_str)
    # If the number is smaller than the position of the missing aminoacid,
    # return the number
    if num <= shift[0]:
        return num_str
    # If the number is bigger than the position of the missing aminoacid,
    # return the number plus the shift
    return str(shift[1] + num)


def ctd_further_check(g, gg):
    return gg['CDS'].qualifiers['systematic_id'][0] in ['SPBC28F2.12'] 


def ctd_convert_to_normal_variant(ctd_substring: str):

    mutations = re.findall(ctd_mutation_regex, ctd_substring)
    deletions= re.findall(ctd_deletion_regex, ctd_substring)
    starting_position = rpb1_dictionary['start']
    repeat_length = len(rpb1_dictionary['positions'])
    out_list = []
    deleted_repeats = list()
    for deletion in deletions:
        # Entire deletion, always correct, and takes over everything else
        if '(' not in deletion:
            return '{}-{}'.format(rpb1_dictionary['start'] + 1, rpb1_dictionary['start'] + rpb1_dictionary['length'])
        # Deletion with repeat number
        match = re.search(ctd_repetition_regex_with_groups, deletion)
        start, stop, step = match.groups()
        step = int(step) if step is not None else 1
        start, stop = int(start), int(stop)
        deleted_ranges = list()

        for repeat_number in range(start, stop + 1, step):
            deleted_repeats.append(repeat_number)
            repeat_start = (repeat_number - 1) * repeat_length + starting_position + 1
            repeat_end = repeat_start + repeat_length - 1
            deleted_ranges.append((repeat_start, repeat_end))
        if step == 1:
            out_list.append('{}-{}'.format(deleted_ranges[0][0], deleted_ranges[-1][1]))
        else:
            out_list.append(','.join(['{}-{}'.format(*x) for x in deleted_ranges]))

    for mutation in mutations:
        original_residue, index_in_repeat, replaced_by = re.search(r'([A-Za-z])(\d+)([A-Za-z])', mutation).groups()
        index_in_repeat = int(index_in_repeat)
        repeats_where_residue_is_present = rpb1_dictionary['positions'][mutation[:2]]
        match = re.search(ctd_repetition_regex_with_groups, mutation)
        if match is None:
            start, stop, step = 1, rpb1_dictionary['nb_repeats'], 1
        else:
            start, stop, step = match.groups()
            step = int(step) if step is not None else 1
            start, stop = int(start), int(stop)
        for repeat_number in range(start, stop + 1, step):
            if repeat_number in deleted_repeats:
                continue
            if repeat_number - 1 not in repeats_where_residue_is_present:
                continue
            mutated_position = index_in_repeat + (repeat_number - 1) * repeat_length + starting_position
            out_list.append('{}{}{}'.format(original_residue, mutated_position, replaced_by))
    out_str = ','.join(out_list)

    return re.sub(r'(\d+)', lambda x: apply_shift(x, rpb1_dictionary['shift']), out_str)


def ctd_check_sequence(ctd_substring: str):
    sequence_errors = []
    mutations = re.findall(ctd_mutation_regex, ctd_substring)

    for mutation in mutations:
        if mutation[:2] not in rpb1_dictionary['positions']:
            sequence_errors.append('CTD-' + mutation[:2])

    match = re.search(ctd_repetition_regex_with_groups, ctd_substring)
    start, stop, step = match.groups()
    if int(stop) > rpb1_dictionary['nb_repeats']:
        sequence_errors.append('CTD-r' + stop)
    if int(start) > rpb1_dictionary['nb_repeats']:
        sequence_errors.append('CTD-r' + start)

    return '/'.join(sequence_errors)


def ctd_format_for_transvar(capture_groups: list[str], gene: dict) -> list[str]:
    ctd_substring = capture_groups[0]
    result = list()
    for ele in ctd_convert_to_normal_variant(ctd_substring).split(','):
        if '-' in ele:
            result.append('p.{}_{}del'.format(*ele.split('-')))
        else:
            result.append('p.{}'.format(ele))
    return result


def ctd_apply_syntax(ctd_substring: str):
    ctd_bits = re.findall(f'({ctd_mutation_regex}|{ctd_deletion_regex})', ctd_substring)
    return 'CTD-' + ','.join(ctd_bits)
