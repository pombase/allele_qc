aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

allowed_types = {
    frozenset({'amino_acid_mutation'}): 'amino_acid_mutation',
    frozenset({'partial_amino_acid_deletion'}): 'partial_amino_acid_deletion',
    frozenset({'amino_acid_mutation', 'partial_amino_acid_deletion'}): 'amino_acid_deletion_and_mutation',
    frozenset({'amino_acid_insertion'}): 'amino_acid_insertion',
    frozenset({'amino_acid_insertion', 'partial_amino_acid_deletion'}): 'amino_acid_insertion_and_deletion',
    frozenset({'amino_acid_insertion', 'amino_acid_mutation'}): 'amino_acid_insertion_and_mutation',
    frozenset({'disruption'}): 'disruption',
    frozenset({'nonsense_mutation'}): 'nonsense_mutation',
    frozenset({'amino_acid_mutation', 'nonsense_mutation'}): 'other',
    frozenset({'unknown'}): 'unknown',
}


def check_position_exists(aa_pos, gene):
    peptide_seq = gene['translation']
    if aa_pos > len(peptide_seq):
        return f'position {aa_pos} does not exist, peptide length is {len(peptide_seq)}'
    return ''


def check_aminoacid_at_pos(aa, aa_pos, gene):

    # Check if the position is valid
    check_pos = check_position_exists(aa_pos, gene)
    if check_pos:
        return check_pos

    # Check if the aminoacid in that position is correct
    # 1-based index
    zero_based_pos = aa_pos - 1
    peptide_seq = gene['translation']
    if peptide_seq[zero_based_pos] == aa:
        return ''

    out_str = f'no {aa} at position {aa_pos}'
    if zero_based_pos + 1 < len(peptide_seq) and peptide_seq[zero_based_pos + 1] == aa:
        out_str += f', but found at {aa_pos + 1}'
    if peptide_seq[zero_based_pos - 1] == aa:
        out_str += f', but found at {aa_pos - 1}'

    return out_str


def check_sequence_single_aa(groups, gene):
    aa = groups[0]
    aa_pos = int(groups[1])
    return check_aminoacid_at_pos(aa, aa_pos, gene)


def check_sequence_multiple_aa(groups, gene):

    pos_first_aa = int(groups[1])
    results_list = list()
    for i, aa in enumerate(groups[0]):
        results_list.append(check_aminoacid_at_pos(aa, pos_first_aa + i, gene))

    output = '|'.join([r for r in results_list if r])
    if len(output):
        return output
    else:
        return ''


def check_multiple_positions(groups, gene):

    results_list = list()
    for pos in groups:
        results_list.append(check_position_exists(int(pos), gene))

    output = '|'.join([r for r in results_list if r])
    if len(output):
        return output
    else:
        return ''


grammar = [
    {
        'type': 'amino_acid_mutation',
        'rule_name': 'single_aa',
        'regex': f'(?<!{aa})({aa})(\d+)({aa})(?!{aa})',
        'apply_syntax': lambda g: ''.join(g).upper(),
        'check_invalid': lambda g: '',
        'check_sequence': check_sequence_single_aa
    },
    {
        'type': 'amino_acid_mutation',
        'rule_name': 'multiple_aa',
        # This is only valid for cases with two aminoacids or more (not to clash with amino_acid_insertion:usual)
        'regex': f'({aa}{aa}+)-?(\d+)-?({aa}+)(?!\d)',
        # We fix the case in which dashes are used for a single aa substitution: K-90-R
        'apply_syntax': lambda g: '-'.join(g).upper() if len(g[0]) != 1 else ''.join(g).upper(),
        'check_invalid': lambda g: f'lengths don\'t match: {g[0]}-{g[2]}' if len(g[0]) != len(g[2]) else '',
        'check_sequence': check_sequence_multiple_aa
    },
    {
        'type': 'nonsense_mutation',
        'rule_name': 'stop_codon_text',
        'regex': f'({aa})(\d+)[^a-zA-Z0-9]*(?i:ochre|stop|amber|opal)',
        'apply_syntax': lambda g: ''.join(g).upper() + '*',
        'check_invalid': lambda g: '',
        'check_sequence': check_sequence_single_aa
    },
    {
        'type': 'nonsense_mutation',
        'rule_name': 'stop_codon_star',
        'regex': f'({aa})(\d+)(\*)',
        'apply_syntax': lambda g: ''.join(g[:2]).upper() + '*',
        'check_invalid': lambda g: '',
        'check_sequence': check_sequence_single_aa
    },
    # {
    #     'type': 'nonsense_mutation',
    #     'rule_name': 'stop_codon_aa_missing',
    #     'regex': f'({aa})(\d+)(\*)',
    #     'apply_syntax': lambda g: ''.join(g[:2]).upper()+'*',
    #     'check_invalid': lambda g: ''
    # },
    {
        'type': 'partial_amino_acid_deletion',
        'rule_name': 'multiple_aa',
        'regex': f'(?<!{aa})(\d+)\s*[-–]\s*(\d+)(?!{aa})(\s+Δaa)?',
        'apply_syntax': lambda g: '-'.join(g[:2]).upper(),
        'check_invalid': lambda g: '',
        'check_sequence': lambda groups, gene: check_multiple_positions(groups[:2], gene)
    },
    {
        'type': 'partial_amino_acid_deletion',
        'rule_name': 'single_aa',
        'regex': f'(?<!{aa})(\d+)(?!{aa})(\s+Δaa)?',
        'apply_syntax': lambda g: g[0],
        'check_invalid': lambda g: '',
        'check_sequence': lambda groups, gene: check_multiple_positions(groups[:1], gene)
    },
    {
        'type': 'amino_acid_insertion',
        'rule_name': 'usual',
        'regex': f'({aa}?)(\d+)-?({aa}+)(?!\d)',
        'apply_syntax': lambda g: '-'.join(g[1:]).upper(),
        'check_invalid': lambda g: '',
        'check_sequence': lambda groups, gene: check_multiple_positions(groups[1:2], gene) if not groups[0] else check_sequence_single_aa(groups, gene)
    },
    {
        'type': 'unknown',
        'rule_name': 'empty',
        'regex': '^$',
        'apply_syntax': lambda g: 'unknown',
        'check_invalid': lambda g: '',
        'check_sequence': lambda g, gg: ''
    },
    {
        'type': 'unknown',
        'rule_name': 'unknown',
        'regex': '^unknown$',
        'apply_syntax': lambda g: 'unknown',
        'check_invalid': lambda g: '',
        'check_sequence': lambda g, gg: ''
    }
]
