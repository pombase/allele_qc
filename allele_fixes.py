import re
from genome_functions import get_other_index_from_alignment
from Bio.Seq import Seq
import regex


def shift_coordinates_by_x(input_str, shift_value):
    """
    Add shift_value to all numbers in input_str.
    """
    splitted_string = re.split(r'(\d+)', input_str)
    splitted_string_new_coords = list()
    for ele in splitted_string:
        # Shift digits only
        if not ele.isdigit():
            splitted_string_new_coords.append(ele)
        else:
            splitted_string_new_coords.append(str(int(ele) + shift_value))
    return ''.join(splitted_string_new_coords)


def position_or_index_exists(target, sequence):
    """
    Check if target (e.g. 'A123', 'P124V', '12-14', '15') exist in sequence. For '12-14', '15' we just check that the coordinates exist.
    """
    re_match = re.match('^([a-zA-Z])(\d+)([a-zA-Z*]?)$', target)
    if re_match is not None:
        # We are dealing with a position (A123) or mutation (A123P)
        # We check that the indicated aminoacid is at that position.
        index = int(re_match.groups()[1]) - 1
        return index < len(sequence) and sequence[index] == target[0]
    else:
        # We are dealing with indexes only, so we just check whether all integers are smaller than the sequence length
        indexes = [(int(i) - 1) for i in re.findall(r'\d+', target)]
        return all(i < len(sequence) for i in indexes)


def change_to_new_indexes_from_old_coordinates(old_alignment, new_alignment, target):
    """
    Substitute all numbers in the string `target` applying get_other_index_from_alignment to each of them.
    If any of the calls to get_other_index_from_alignment returns None, this function returns None.
    """

    # Split by number, keeping the number in (using the capture group)
    splitted_string = re.split(r'(\d+)', target)
    splitted_string_new_coords = list()
    for ele in splitted_string:
        # Apply get_other_index_from_alignment to the digits only
        if not ele.isdigit():
            splitted_string_new_coords.append(ele)
        else:
            new_value_zero_based = get_other_index_from_alignment(old_alignment, new_alignment, int(ele) - 1)
            if new_value_zero_based is None:
                return None
            splitted_string_new_coords.append(str(new_value_zero_based + 1))
    return ''.join(splitted_string_new_coords)


def old_coords_fix(coordinate_changes, targets):
    """
    Propose a fix if the sequence positions and coordinates proposed in targets (e.g. ['A123', 'P124V', '12-14'])
    match those of an old sequence (for '12-14' we just check that the coordinate exists), and for the rest we
    check that the indicated aminoacid is there in the old and new sequence. (see test_coordinate_change)
    """
    out_list = list()
    for prev_coord in coordinate_changes:

        new_alignment = prev_coord['new_alignment']
        old_alignment = prev_coord['old_alignment']
        new_sequence = new_alignment.replace('-', '')
        old_sequence = old_alignment.replace('-', '')

        this_revision = {'revision': prev_coord['revision'], 'location': prev_coord['old_coord'], 'values': list()}
        # remap the coordinates
        for t in targets:
            # The position must exist in the old sequence
            if not position_or_index_exists(t, old_sequence):
                break

            # We now remap from old coordinates to new ones. It may return None if the position in the old sequence
            # does not exist in the new one. We also checked whether the new value matches the new sequence.
            target_remapped = change_to_new_indexes_from_old_coordinates(old_alignment, new_alignment, t)
            if target_remapped is None or not position_or_index_exists(target_remapped, new_sequence):
                break

            this_revision['values'].append(target_remapped)
        else:
            # Only append the value if break was not triggered
            this_revision['values'] = ','.join(this_revision['values'])
            out_list.append(this_revision)

    return out_list


def multi_shift_fix(seq, targets):
    """
    Check if changing all coordinates by a fixed amount they match the sequence. Return
    all possible solutions

    seq = AVPPAVPPP
    targets = ['A3', 'V4'] # or ['A3L', 'V4L']
    returns ['A1,V2', 'A5,V6']

    It also works with only index targets, checking if the sequence is shorter than the
    indicated value.
    seq = AVPPAVPPP
    targets = ['A3', 'V4', '8']
    returns ['A1,V2'] # Because position 10 does not exist
    """

    all_indexes = list()
    for target in targets:
        all_indexes.extend([(int(i) - 1) for i in re.findall(r'\d+', target)])

    minus_shift = - min(all_indexes)
    plus_shift = len(seq) + minus_shift
    out_list = list()

    for shift_amount in range(minus_shift, plus_shift):
        shifted_targets = [shift_coordinates_by_x(target, shift_amount) for target in targets]
        if all(position_or_index_exists(shifted_target, seq) for shifted_target in shifted_targets):
            out_list.append(','.join(shifted_targets))

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
