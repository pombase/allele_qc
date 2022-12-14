from genome_functions import get_nt_at_genome_position, gene_coords2genome_coords


allowed_types = {
    frozenset({'amino_acid_mutation'}): 'amino_acid_mutation',
    frozenset({'partial_amino_acid_deletion'}): 'partial_amino_acid_deletion',
    frozenset({'amino_acid_mutation', 'partial_amino_acid_deletion'}): 'amino_acid_deletion_and_mutation',
    frozenset({'amino_acid_insertion'}): 'amino_acid_insertion',
    frozenset({'amino_acid_insertion', 'partial_amino_acid_deletion'}): 'amino_acid_insertion_and_deletion',
    frozenset({'amino_acid_insertion', 'amino_acid_mutation'}): 'amino_acid_insertion_and_mutation',
    frozenset({'disruption'}): 'disruption',
    frozenset({'nonsense_mutation'}): 'nonsense_mutation',
    frozenset({'amino_acid_mutation', 'nonsense_mutation'}): 'amino_acid_other',
    frozenset({'nucleotide_mutation'}): 'nucleotide_mutation',
    frozenset({'nucleotide_insertion'}): 'nucleotide_insertion',
    frozenset({'partial_nucleotide_deletion'}): 'partial_nucleotide_deletion',
    frozenset({'nonsense_mutation', 'amino_acid_insertion'}): 'amino_acid_other'
}


def check_position_doesnt_exist(pos, gene, seq_type):
    """
    Return error string if the position is beyond the end of the sequence.
    """
    # position zero is not allowed for 1-based indexing of peptides
    if pos == 0:
        return 'position 0 is not allowed when 1-based indexing'
    if seq_type == 'peptide':
        peptide_seq = gene['peptide']
        if pos > len(peptide_seq):
            return f'position {pos} does not exist, peptide length is {len(peptide_seq)}'
        return ''

    try:
        gene_coords2genome_coords(pos, gene)
    except IndexError:
        return f'cannot access genome position {pos}'
    except ValueError as e:
        return e.args[0]


def check_value_at_pos(value, pos, gene, seq_type, append_suggestion=True):
    """
    Return error string if the position is beyond the end of the sequence or the indicated
    value is not at that position.
    """
    # Check if the position is valid
    check_pos = check_position_doesnt_exist(pos, gene, seq_type)
    if check_pos:
        return check_pos
    if seq_type == 'dna':
        value = value.replace('u', 't').replace('U', 'T')
    # Check if the value in that position is correct
    if seq_type == 'peptide':
        def get_value_at_pos(p):
            return gene['peptide'][p - 1]
    else:
        def get_value_at_pos(p):
            return get_nt_at_genome_position(p, gene, gene['contig'])

    if get_value_at_pos(pos) == value:
        return ''

    out_str = f'{value}{pos}'

    if append_suggestion:
        for i in [1, -1]:
            if not check_position_doesnt_exist(pos + i, gene, seq_type) and (get_value_at_pos(pos + i) == value):
                out_str += f'>{value}{pos + i}'

    return out_str


def check_sequence_single_pos(groups, gene, seq_type, append_suggestion=True):
    """
    Check that a single position in the sequence exists, defined in regex capture groups. In `groups` the first
    element is the aminoacid / nucleotide and the second is the number, e.g. in V123A, groups are ['V', '123','A']
    """
    value = groups[0]
    pos = int(groups[1])
    return check_value_at_pos(value, pos, gene, seq_type, append_suggestion)


def check_sequence_multiple_pos(groups, gene, seq_type, append_suggestion=True):
    """
    Check multiple positions in the sequence from the capture groups of the regex.
    This is called for example in 'VPL-234-AAA'.
    The capture groups from the grammar VPL and 234 in the above example. If several errors are found, they are returned as `\` separated
    strings. E.g. `V234/P235/L236`.
    """

    pos_first = int(groups[1])
    results_list = list()
    # Iterate over chars of string
    for i, value in enumerate(groups[0]):
        results_list.append(check_value_at_pos(value, pos_first + i, gene, seq_type, append_suggestion))

    output = '/'.join([r for r in results_list if r])
    if len(output):
        return output
    else:
        return ''


def check_multiple_positions_dont_exist(groups, gene, seq_type):

    results_list = list()
    for pos in groups:
        results_list.append(check_position_doesnt_exist(int(pos), gene, seq_type))

    output = '/'.join([r for r in results_list if r])
    if len(output):
        return output
    else:
        return ''


# Variable with all aminoacids, to be used in the regex
aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

aminoacid_grammar = [
    {
        'type': 'amino_acid_mutation',
        'rule_name': 'single_aa',
        'regex': f'(?<!{aa})({aa})(\d+)({aa})(?!{aa})',
        'apply_syntax': lambda g: ''.join(g).upper(),
        'check_invalid': lambda g: '',
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'peptide'),
        'coordinate_indexes': (1,)
    },
    {
        'type': 'amino_acid_mutation',
        'rule_name': 'multiple_aa',
        # This is only valid for cases with two aminoacids or more (not to clash with amino_acid_insertion:usual)
        'regex': f'(?<!\d)({aa}{aa}+)-?(\d+)-?({aa}+)(?!\d)',
        'apply_syntax': lambda g: '-'.join(g).upper(),
        'check_invalid': lambda g: f'lengths don\'t match: {g[0]}-{g[2]}' if len(g[0]) != len(g[2]) else '',
        'check_sequence': lambda g, gg: check_sequence_multiple_pos(g, gg, 'peptide'),
        'coordinate_indexes': (1,)
    },
    {
        'type': 'nonsense_mutation',
        'rule_name': 'stop_codon_text',
        'regex': f'({aa})(\d+)[^a-zA-Z0-9]*(?i:ochre|stop|amber|opal)',
        'apply_syntax': lambda g: ''.join(g).upper() + '*',
        'check_invalid': lambda g: '',
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'peptide'),
        'coordinate_indexes': (1,)
    },
    {
        'type': 'nonsense_mutation',
        'rule_name': 'stop_codon_star',
        'regex': f'({aa})(\d+)(\*)',
        'apply_syntax': lambda g: ''.join(g[:2]).upper() + '*',
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'peptide'),
        'coordinate_indexes': (1,)
    },
    {
        'type': 'partial_amino_acid_deletion',
        'rule_name': 'multiple_aa',
        'regex': f'(?<!{aa})(\d+)\s*[-–]\s*(\d+)(?!{aa})(\s+Δaa)?',
        'apply_syntax': lambda g: '-'.join(sorted(g[:2], key=int)).upper(),
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups[:2], gene, 'peptide'),
        'coordinate_indexes': (0, 1)
    },
    {
        'type': 'partial_amino_acid_deletion',
        'rule_name': 'single_aa',
        'regex': f'(?<!{aa})(\d+)(?!{aa})(\s+Δaa)?',
        'apply_syntax': lambda g: g[0],
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups[:1], gene, 'peptide'),
        'coordinate_indexes': (0,)
    },
    {
        'type': 'amino_acid_insertion',
        'rule_name': 'usual',
        'regex': f'({aa}?)(\d+)-?({aa}+)(?!\d)',
        'apply_syntax': lambda g: '-'.join(g[1:]).upper(),
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups[1:2], gene, 'peptide') if not groups[0] else check_sequence_single_pos(groups, gene, 'peptide'),
        'coordinate_indexes': (1,)
    }
]


def format_negatives(input_list, indexes):
    output_list = list(input_list[:])
    for index in indexes:
        output_list[index] = output_list[index] if int(output_list[index]) > 0 else f'({output_list[index]})'
    return output_list


# Variable with all the nucleotides, to be used in the regex
# We allow the U, but we will replace it by T
nt = 'ACGUT'
nt = nt + nt.lower()
nt = f'[{nt}]'

nucleotide_grammar = [
    {
        'type': 'nucleotide_mutation',
        'rule_name': 'single_nt',
        # Negative numbers are common
        'regex': f'(?<!{nt})({nt})(-?\d+)({nt})(?!{nt})',
        'apply_syntax': lambda g: ''.join(format_negatives(g, [1])).upper().replace('U', 'T'),
        'check_invalid': lambda g: '',
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'dna')
    },
    {
        'type': 'nucleotide_mutation',
        'rule_name': 'multiple_nt',
        # This is only valid for cases with two nts or more (not to clash with nucleotide_insertion:usual)
        # Note the non-greedy flanking dashes, to prioritise the dash for negative numbers
        'regex': f'({nt}{nt}+)-??(-?\d+)-??({nt}+)(?!\d)',
        'apply_syntax': lambda g: ('-'.join(format_negatives(g, [1])) if len(g[0]) != 1 else ''.join(g)).upper().replace('U', 'T'),
        'check_invalid': lambda g: f'lengths don\'t match: {g[0]}-{g[2]}' if len(g[0]) != len(g[2]) else '',
        'check_sequence': lambda g, gg: check_sequence_multiple_pos(g, gg, 'dna')
    },
    {
        'type': 'partial_nucleotide_deletion',
        'rule_name': 'usual',
        'regex': f'(?<!{nt})(-?\d+)\s*[-–]\s*(-?\d+)(?!{nt})',
        'apply_syntax': lambda g: '-'.join(format_negatives(sorted(g, key=int), [0, 1])).upper(),
        'check_invalid': lambda g: '',
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups, gene, 'dna')
    },

]

disruption_grammar = [
    {
        'type': 'disruption',
        'rule_name': 'usual',
        'regex': '([a-zA-Z]{3}\d+|SP[A-Z0-9]+\.[A-Za-z0-9]+)::(.+?)\+?\s*$',
        'apply_syntax': lambda g: f'{g[0]}::{g[1]}',
    }
]
