import re
from genome_functions import get_other_index_from_alignment
from Bio.Seq import Seq
import regex


# def read_sequence_position_from_string(input_str: str):


def extract_groups_from_targets(targets: list[str]):
    """
    From targets ['V123A', 'P125L']
    Returns list of dicts:
    [
        {'value': 'V', 'number': 122, 'new_value': 'A', 'number_zero': 0},
        {'value': 'P', 'number': 124, 'new_value': 'L', 'number_zero': 2}
    ]
    It sorts by number, and number_zero is substracting result[0]['number'] from all numbers
    Also, note that it uses zero-based indexing in number
    """
    groups = list()
    for t in targets:
        this_dict = dict()
        this_dict['value'], this_dict['number'], this_dict['new_value'] = re.match('^([a-zA-Z])(\d+)([a-zA-Z*]?)$', t).groups()
        this_dict['number'] = int(this_dict['number']) - 1
        groups.append(this_dict)

    groups.sort(key=lambda x: x['number'])

    for group in groups:
        group['number_zero'] = group['number'] - groups[0]['number']

    return groups


def shift_coordinates_by_x(groups, seq, number_key, i):
    """
    Shift the coordinates in groups (see extract_groups_from_targets) by a certain
    index i. `number_key` can be indicated to shift from `number_zero` or `number`.

    Returns the shifted coordinates in new coords, and whether they match the specified
    sequence (seq argument).

    """
    shifted_coords_match = list()
    new_coords = list()
    for group in groups:
        n, value, new_value = (group[number_key], group['value'], group['new_value'],)
        shifted_coords_match.append(value == seq[n + i])
        new_coords.append(value + str(n + i + 1) + new_value)

    return new_coords, shifted_coords_match


def old_coords_fix(coordinate_changes, targets):
    out_list = list()
    for prev_coord in coordinate_changes:

        new_alignment = prev_coord['new_alignment']
        old_alignment = prev_coord['old_alignment']
        revision = prev_coord['revision']
        current_seq = new_alignment.replace('-', '')

        numbers = list()
        values = list()
        for t in targets:
            numbers.append(int(re.search(r'\d+', t).group()) - 1)
            values.append(t[0])

        this_revision = {'revision': revision, 'location': prev_coord['old_coord'], 'values': list(), 'matches': list()}

        for v, old_index in zip(values, numbers):
            new_index = get_other_index_from_alignment(old_alignment, new_alignment, old_index)
            if new_index is None:
                break
            this_revision['values'].append(f'{current_seq[new_index]}{new_index+1}')
            this_revision['matches'].append(v == current_seq[new_index])
        else:
            # Only append the value if break was not triggered (if any index is None)
            out_list.append(this_revision)

    return out_list


def multi_shift_fix(seq, targets, shift_amount=None):

    groups = extract_groups_from_targets(targets)

    out_list = list()
    for i in range(len(seq) - groups[-1]['number_zero']):
        new_coords, shifted_coords_match = shift_coordinates_by_x(groups, seq, 'number_zero', i)
        if all(shifted_coords_match):
            out_list.append(','.join(new_coords))

    return out_list


def primer_mutagenesis(main_seq, primer_seq, allowed_mismatches, has_peptide):
    main_seq = str(main_seq)
    primer_seq = Seq(primer_seq)
    out_str = ''

    for reverse_the_primer in [True, False]:
        for i in range(1, allowed_mismatches):
            primer = str(primer_seq) if not reverse_the_primer else str(primer_seq.reverse_complement())
            pattern = f'({primer})' + '{s<=' + str(i) + '}'
            matches = list(regex.finditer(pattern, main_seq))

            for match in matches:
                match: regex.Match
                start = match.start()
                end = start + len(primer)
                primer_DNA = main_seq[:start] + primer + main_seq[end:]
            if len(matches) > 0:
                break
        else:
            continue

        matches_rvs = ''
        mutations = ''
        for i in range(len(main_seq)):
            matches_rvs = matches_rvs + (' ' if primer_DNA[i] == main_seq[i] else '|')
            if primer_DNA[i] != main_seq[i]:
                mutations += f'{main_seq[i]}{i+1}{primer_DNA[i]},'

        if reverse_the_primer:
            out_str = out_str + f'> with the reverse primer and {i} mistmatches\n'
        else:
            out_str = out_str + f'> with the forward primer and {i} mistmatches\n'

        out_str = out_str + f'mutations observed at DNA level: {mutations}\n'
        out_str = out_str + f'original_sequence:  {main_seq}\n'
        out_str = out_str + f'                    {matches_rvs}\n'
        out_str = out_str + f'mutated_sequence:   {primer_DNA}\n'

        if not has_peptide:
            continue

        peptide_seq = Seq(main_seq).translate()
        peptide_rvs = Seq(primer_DNA).translate()

        matches_rvs = ''
        mutations = ''
        for i in range(len(peptide_seq)):
            matches_rvs = matches_rvs + \
                (' ' if peptide_rvs[i] == peptide_seq[i] else '|')
            if peptide_rvs[i] != peptide_seq[i]:
                mutations += f'{peptide_seq[i]}{i+1}{peptide_rvs[i]},'

        out_str = out_str + '\n'
        out_str = out_str + f'mutations observed at peptide level:{mutations}\n'
        out_str = out_str + f'original_sequence: {peptide_seq}\n'
        out_str = out_str + f'                   {matches_rvs}\n'
        out_str = out_str + f'mutated_sequence:  {peptide_rvs}\n'

    if not out_str:
        return 'Nothing found :_(. This works with cDNA, perhaps primer aligns at exon border?'

    return out_str
