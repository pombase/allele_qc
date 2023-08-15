import re

aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

ctd_repetition_regex = r'\(r\d+-r\d+(?:-\d+)?\)'
ctd_repetition_regex_with_groups = r'\(r(\d+)-r(\d+)(?:-(\d+))?\)'
ctd_mutation_regex = f'{aa}\d+{aa}(?:{ctd_repetition_regex})?'
ctd_deletion_regex = f'(?:delta|\u0394)(?:{ctd_repetition_regex})?'

# This dictionary has keys that are the systematic IDs of the genes that have CTDs,
# the value is a dictionary with:
# - positions: a dictionary with keys that are the residues of the canonical repeat
#   e.g. YSPTSPS for SPBC28F2.12, and the value is the repeats in which the residue is
#   present.
# - shift: a list of integers in which each value represents the shift of residues
#   associated at a given repeat. For instance, see the value -1 set after repeat 3 in SPBC28F2.12
#   that has one amino acid less than the other repeats.
# - start: the position of the first residue of the CTD in the protein (zero-based index)
# - length: the length of the CTD

ctd_dictionary = {
    'SPBC28F2.12':
    {
        'positions': {
            'Y1': range(29),
            'S2': [1, 2] + list(range(4, 29)),
            'P3': [1, ] + list(range(3, 8)) + list(range(9, 29)),
            'T4': [0, 2] + list(range(4, 29)),
            'S5': range(29),
            'P6': range(29),
            'S7': [0, 3] + list(range(4, 29)),
        },
        'shift': [0, 0] + [-1] * 27,
        'start': 1550,
        'length': 202,
    }
}


def ctd_further_check(g, gg):
    return gg['CDS'].qualifiers['systematic_id'][0] in ['SPAC23C4.19', 'SPBC28F2.12'] 


def ctd_convert_to_normal_variant(ctd_substring: str, systematic_id: str):
    mutations = re.findall(ctd_mutation_regex, ctd_substring)
    deletions= re.findall(ctd_deletion_regex, ctd_substring)
    ctd_dict = ctd_dictionary[systematic_id]
    out_list = []
    for deletion in deletions:
        # Entire deletion, always correct, and takes over everything else
        if '(' not in deletion:
            return '{}-{}'.format(ctd_dict['start'] + 1, ctd_dict['start'] + ctd_dict['length'])
        # Deletion with repeat number
        match = re.search(ctd_repetition_regex_with_groups, deletion)
        start, stop, step = match.groups()
        step = int(step) if step is not None else 1
        start, stop = int(start), int(stop)
        deleted_ranges = list()
        starting_position = ctd_dict['start']
        repeat_length = len(ctd_dict['positions'])
        shift = ctd_dict['shift']
        for repeat_number in range(start, stop + 1, step):
            repeat_start = (repeat_number - 1) * repeat_length + starting_position + shift[repeat_number - 1] + 1
            repeat_end = repeat_start + repeat_length - 1
            deleted_ranges.append((repeat_start, repeat_end))
        if step == 1:
            out_list.append('{}-{}'.format(deleted_ranges[0][0], deleted_ranges[-1][1]))
        else:
            out_list.append(','.join(['{}-{}'.format(*x) for x in deleted_ranges]))

    return ','.join(out_list)

def ctd_check_sequence(ctd_substring: str, systematic_id: str):
    sequence_errors = []
    mutations = re.findall(ctd_mutation_regex, ctd_substring)
    deletions = re.findall(ctd_deletion_regex, ctd_substring)

    for deletion in deletions:
        # Entire deletion, always correct
        if '(' not in deletion:
            continue
        # Deletion with repeat number
        match = re.search(ctd_repetition_regex_with_groups)
        start, stop, step = match.groups()
        if int(stop) > len(ctd_dictionary[systematic_id]['shift']):
            sequence_errors.append('CTD-r' + stop)
    for mutation in mutations:
        if mutation[:2] not in ctd_dictionary[systematic_id]['positions']:
            sequence_errors.append('CTD-' + mutation[:2])

    return '/'.join(sequence_errors)



def ctd_apply_syntax(ctd_substring: str):
    ctd_bits = re.findall(f'({ctd_mutation_regex}|{ctd_deletion_regex})', ctd_substring)
    return 'CTD-' + ','.join(ctd_bits)