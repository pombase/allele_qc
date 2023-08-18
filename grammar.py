from genome_functions import get_nt_at_gene_coord, gene_coords2genome_coords
import copy
from Bio.Seq import reverse_complement
from ctd_support import ctd_apply_syntax, ctd_convert_to_normal_variant, ctd_further_check, ctd_format_for_transvar
from ctd_support import ctd_mutation_regex, ctd_deletion_regex

allowed_types_dict = {
    frozenset({'amino_acid_mutation'}): 'amino_acid_mutation',
    frozenset({'partial_amino_acid_deletion'}): 'partial_amino_acid_deletion',
    frozenset({'amino_acid_mutation', 'partial_amino_acid_deletion'}): 'amino_acid_deletion_and_mutation',
    frozenset({'amino_acid_insertion'}): 'amino_acid_insertion',
    frozenset({'amino_acid_insertion', 'partial_amino_acid_deletion'}): 'amino_acid_insertion_and_deletion',
    frozenset({'amino_acid_insertion', 'amino_acid_mutation'}): 'amino_acid_insertion_and_mutation',
    frozenset({'disruption'}): 'disruption',
    frozenset({'nonsense_mutation'}): 'partial_amino_acid_deletion',
    frozenset({'amino_acid_mutation', 'nonsense_mutation'}): 'amino_acid_deletion_and_mutation',
    frozenset({'nucleotide_mutation'}): 'nucleotide_mutation',
    frozenset({'nucleotide_insertion'}): 'nucleotide_insertion',
    frozenset({'partial_nucleotide_deletion'}): 'partial_nucleotide_deletion',
    frozenset({'nonsense_mutation', 'amino_acid_insertion'}): 'amino_acid_insertion_and_deletion',
    frozenset({'nucleotide_mutation', 'partial_nucleotide_deletion'}): 'nucleotide_deletion_and_mutation',
    frozenset({'nucleotide_mutation', 'nucleotide_insertion'}): 'nucleotide_insertion_and_mutation',
    frozenset({'nucleotide_insertion', 'partial_nucleotide_deletion', 'nucleotide_mutation'}): 'nucleotide_insertion_and_deletion_and_mutation',
}

composed_types_dict = {
    'amino_acid_insertion_and_mutation': ['amino_acid_insertion', 'amino_acid_mutation'],
    'amino_acid_deletion_and_mutation': ['amino_acid_mutation', 'partial_amino_acid_deletion'],
    'nucleotide_insertion_and_mutation': ['nucleotide_insertion', 'nucleotide_mutation'],
    'nucleotide_deletion_and_mutation': ['nucleotide_mutation', 'partial_nucleotide_deletion'],
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
            return str(pos)
        return ''

    try:
        gene_coords2genome_coords(pos, gene)
    except IndexError:
        return f'cannot access genome position {pos}'
    except ValueError as e:
        return e.args[0]


def check_value_at_pos(indicated_value, pos, gene, seq_type):
    """
    Return error string if the position is beyond the end of the sequence or the indicated
    indicated_value is not at that position.
    """
    # Check if the position is valid
    check_pos = check_position_doesnt_exist(pos, gene, seq_type)
    if check_pos:
        if pos < 0:
            return f'{indicated_value}({pos})'
        return f'{indicated_value}{pos}'
    if seq_type == 'dna':
        indicated_value = indicated_value.replace('u', 't').replace('U', 'T')
    # Check if the value in that position is correct
    if seq_type == 'peptide':
        def get_value_at_pos(p):
            return gene['peptide'][p - 1]
    else:
        def get_value_at_pos(p):
            return get_nt_at_gene_coord(p, gene, gene['contig'])

    # Very important to be case-insensitive
    if get_value_at_pos(pos).upper() == indicated_value.upper():
        return ''

    if pos < 0:
        return f'{indicated_value}({pos})'
    return f'{indicated_value}{pos}'


def check_sequence_single_pos(groups: list[str], gene, seq_type):
    """
    Check that a single position in the sequence exists, defined in regex capture groups. In `groups` the first
    element is the aminoacid / nucleotide and the second is the number, e.g. in V123A, groups are ['V', '123','A']
    """
    value = groups[0]
    pos = int(groups[1].replace('(', '').replace(')', ''))
    return check_value_at_pos(value, pos, gene, seq_type)


def check_sequence_multiple_pos(groups, gene, seq_type):
    """
    Check multiple positions in the sequence from the capture groups of the regex.
    This is called for example in 'VPL-234-AAA'.
    The capture groups from the grammar VPL and 234 in the above example. If several errors are found, they are returned as `\` separated
    strings. E.g. `V234/P235/L236`.
    """

    pos_first = int(groups[1].replace('(', '').replace(')', ''))
    results_list = list()
    # Iterate over chars of string
    for i, value in enumerate(groups[0]):
        results_list.append(check_value_at_pos(value, pos_first + i, gene, seq_type))

    output = '/'.join([r for r in results_list if r])
    if len(output):
        return output
    else:
        return ''


def check_multiple_positions_dont_exist(groups, gene, seq_type):

    results_list = list()
    for pos in groups:
        results_list.append(check_position_doesnt_exist(int(pos.replace('(', '').replace(')', '')), gene, seq_type))

    output = '/'.join([r for r in results_list if r])
    if len(output):
        return output
    else:
        return ''


# Variable with all aminoacids, to be used in the regex
aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

aminoacid_grammar_old = [
    {
        'type': 'amino_acid_mutation',
        'rule_name': 'single_aa',
        'regex': f'(?<=\\b)({aa})(\d+)({aa})(?=\\b)',
        'apply_syntax': lambda g: ''.join(g).upper(),
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'peptide'),
    },
    {
        'type': 'amino_acid_mutation',
        'rule_name': 'multiple_aa',
        # This is only valid for cases with two aminoacids or more (not to clash with amino_acid_insertion)
        'regex': f'(?<=\\b)({aa}{aa}+)-?(\d+)-?({aa}+)(?=\\b)',
        'apply_syntax': lambda g: '-'.join(g).upper(),
        'check_sequence': lambda g, gg: check_sequence_multiple_pos(g, gg, 'peptide'),
    },
    {
        'type': 'nonsense_mutation',
        'rule_name': 'stop_codon_text',
        'regex': f'({aa})(\d+)[^a-zA-Z0-9]*(?i:ochre|stop|amber|opal)',
        'apply_syntax': lambda g: ''.join(g).upper() + '*',
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'peptide'),
    },
    {
        'type': 'nonsense_mutation',
        'rule_name': 'stop_codon_star',
        'regex': f'({aa})(\d+)(\*)',
        'apply_syntax': lambda g: ''.join(g[:2]).upper() + '*',
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'peptide'),
    },
    {
        'type': 'partial_amino_acid_deletion',
        'rule_name': 'multiple_aa',
        'regex': f'(?<!{aa})(\d+)\s*[-–]\s*(\d+)(?!{aa})(?:\s+Δaa)?',
        'apply_syntax': lambda g: '-'.join(sorted(g, key=int)).upper(),
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups, gene, 'peptide'),
    },
    {
        'type': 'partial_amino_acid_deletion',
        'rule_name': 'single_aa',
        'regex': f'(?<!{aa})(\d+)(?!{aa})(?:\s+Δaa)?',
        'apply_syntax': lambda g: g[0],
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups, gene, 'peptide'),
    },
    # We split the insertion into two cases, one where a single aminoacid is inserted, in which the dash
    # is compulsory, and one where the dash is optional, for more than one. Otherwise A123V would match
    # this and the amino_acid_mutation.
    {
        'type': 'amino_acid_insertion',
        'rule_name': 'single',
        'regex': f'({aa})(\d+)-({aa})(?=\\b)',
        'apply_syntax': lambda g: f'{g[0]}{g[1]}-{g[2]}'.upper(),
        'check_sequence': lambda groups, gene: check_sequence_single_pos(groups, gene, 'peptide'),
    },
    {
        'type': 'amino_acid_insertion',
        'rule_name': 'multiple',
        'regex': f'({aa})(\d+)-?({aa}{aa}+)(?=\\b)',
        'apply_syntax': lambda g: f'{g[0]}{g[1]}-{g[2]}'.upper(),
        'check_sequence': lambda groups, gene: check_sequence_single_pos(groups, gene, 'peptide'),
    }
]


def format_negatives(input_list: list[str], indexes: list[int]):
    output_list = list(input_list[:])
    for index in indexes:
        this_value = output_list[index].replace('(', '').replace(')', '')
        output_list[index] = this_value if int(this_value) > 0 else f'({this_value})'
    return output_list


# Variable with all the nucleotides, to be used in the regex
# We allow the U, but we will replace it by T
nt = 'ACGUT'
nt = nt + nt.lower()
nt = f'[{nt}]'

# We favour formatting negative numbers with parenthesis
# this regex captures both positive and negative numbers without parenthesis, and
# negative numbers with parenthesis
num = '(\(-\d+\)|(?<!\()-?\d+(?!\)))'

nucleotide_grammar_old = [
    {
        'type': 'nucleotide_mutation',
        'rule_name': 'single_nt',
        # Negative numbers are common
        'regex': f'(?<=\\b)({nt}){num}({nt})(?=\\b)',
        'apply_syntax': lambda g: ''.join(format_negatives(g, [1])).upper().replace('U', 'T'),
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'dna')
    },
    {
        'type': 'nucleotide_mutation',
        'rule_name': 'multiple_nt',
        # This is only valid for cases with two nts or more (not to clash with nucleotide_insertion:usual)
        # Cases contemplated here:
        # AA-23-TT (positive number correctly formatted)
        # AA-(-23)-TT (negative number correctly formatted)
        # AA23TT (positive number without dashes)
        # AA-23TT (negative number without dashes nor parenthesis)
        # AA(-23)TT (negative number without dashes)
        # AA--23-TT (negative number without parenthesis)
        # Note the use of positive and negative lookahead / lookbehind for dashes to include both cases
        'regex': f'({nt}{nt}+)-?((?<=-)(?:-?\d+|\(-\d+\))(?=-)|(?<!-)(?:-?\d+|\(-\d+\))(?!-))-?({nt}+)(?=\\b)',
        'apply_syntax': lambda g: ('-'.join(format_negatives(g, [1]))).upper().replace('U', 'T'),
        'check_sequence': lambda g, gg: check_sequence_multiple_pos(g, gg, 'dna')
    },
    {
        'type': 'partial_nucleotide_deletion',
        'rule_name': 'usual',
        'regex': f'(?<!{nt}){num}\s*[-–]\s*{num}(?!{nt})',
        'apply_syntax': lambda g: '-'.join(format_negatives(sorted(g, key=lambda x: int(x.replace('(', '').replace(')', ''))), [0, 1])).upper(),
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups, gene, 'dna')
    },
    {
        'type': 'partial_nucleotide_deletion',
        'rule_name': 'single_nt',
        'regex': f'(?<!{nt}){num}(?!{nt})',
        'apply_syntax': lambda g: format_negatives(g, [0])[0],
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups[:1], gene, 'dna'),
    },
    # We split the insertion into two cases, one where a single nt is inserted, in which the dash
    # is compulsory, and one where the dash is optional, for more than one. Otherwise A123T would match
    # this and the nucleotide_mutation.
    {
        'type': 'nucleotide_insertion',
        'rule_name': 'single',
        'regex': f'({nt}){num}-({nt})(?=\\b)',
        'apply_syntax': lambda g: f'{g[0]}{format_negatives(g[1:2],[0])[0]}-{g[2]}'.upper().replace('U', 'T'),
        'check_sequence': lambda groups, gene: check_sequence_single_pos(groups, gene, 'dna'),
    },
    {
        'type': 'nucleotide_insertion',
        'rule_name': 'multiple',
        'regex': f'({nt}){num}-?({nt}{nt}+)(?=\\b)',
        'apply_syntax': lambda g: f'{g[0]}{format_negatives(g[1:2],[0])[0]}-{g[2]}'.upper().replace('U', 'T'),
        'check_sequence': lambda groups, gene: check_sequence_single_pos(groups, gene, 'dna'),
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


# New grammars - here there are a lot of re-used regex, so we use variables to avoid repetition

multi_aa_regex = f'(?<=\\b)({aa}+)-?(\d+)-?({aa}+)(?=\\b)'
multi_aa_apply_syntax = lambda g: ''.join(g).upper()
multi_aa_check_sequence = lambda g, gg: check_sequence_multiple_pos(g, gg, 'peptide')


def multi_aa_format_for_transvar(g: list[str], gg=None) -> list[str]:
    # First aminoacid + first position
    first_residue = g[0][0] + g[1]
    # Last aminoacid + last position
    last_residue = g[0][-1] + str(int(g[1]) + len(g[0]) - 1)

    return [f'p.{first_residue}_{last_residue}delins{g[2]}']


def insertion_aa_format_for_transvar(g: list[str], gg=None) -> list[str]:
    # Last aminoacid + last position
    last_pos = int(g[1]) + len(g[0]) - 1
    last_residue = g[0][-1] + str(last_pos)

    return [f'p.{last_residue}_{last_pos + 1}ins{g[2][1:]}']


aminoacid_grammar = [
    {
        'type': 'amino_acid_mutation',
        'rule_name': 'single_aa',
        'regex': f'(?<=\\b)({aa})(\d+)({aa})(?=\\b)',
        'apply_syntax': lambda g: ''.join(g).upper(),
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'peptide'),
        'further_check': lambda g, gg: g[0] != g[2],
        'format_for_transvar': lambda g, gg: [f'p.{g[0]}{g[1]}{g[2]}']
    },
    {
        'type': 'amino_acid_mutation',
        'rule_name': 'multiple_aa',
        'regex': multi_aa_regex,
        'apply_syntax': multi_aa_apply_syntax,
        'check_sequence': multi_aa_check_sequence,
        # It is only a mutation if the number of aminoacids before and after is the same
        'further_check': lambda g, gg: (len(g[0]) == len(g[2])) & (g[0] != g[2]),
        'format_for_transvar': multi_aa_format_for_transvar
    },
    {
        'type': 'amino_acid_deletion_and_mutation',
        'rule_name': 'multiple_aa',
        'regex': multi_aa_regex,
        'apply_syntax': multi_aa_apply_syntax,
        'check_sequence': multi_aa_check_sequence,
        # TODO: Here we could even check that a partial deletion has not been written using this syntax: e.g. AVTGLA123AA, but probably rare enough.
        'further_check': lambda g, gg: len(g[0]) > len(g[2]),
        'format_for_transvar': multi_aa_format_for_transvar
    },
    {
        'type': 'amino_acid_insertion_and_mutation',
        'rule_name': 'multiple_aa',
        'regex': multi_aa_regex,
        'apply_syntax': multi_aa_apply_syntax,
        'check_sequence': multi_aa_check_sequence,
        # We don't want to account insertions here
        'further_check': lambda g, gg: (len(g[0]) < len(g[2])) & (not g[2].startswith(g[0])),
        'format_for_transvar': multi_aa_format_for_transvar
    },
    {
        'type': 'amino_acid_insertion',
        'rule_name': 'standard',
        'regex': multi_aa_regex,
        'apply_syntax': multi_aa_apply_syntax,
        'check_sequence': multi_aa_check_sequence,
        'further_check': lambda g, gg: (len(g[0]) < len(g[2])) & (g[2].startswith(g[0])),
        'format_for_transvar': insertion_aa_format_for_transvar
    },
    {
        'type': 'nonsense_mutation',
        'rule_name': 'stop_codon_text',
        'regex': f'({aa})(\d+)[^a-zA-Z0-9]*(?i:ochre|stop|amber|opal)',
        'apply_syntax': lambda g: ''.join(g).upper() + '*',
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'peptide'),
        'format_for_transvar': lambda g, gg: [f'p.{g[0]}{g[1]}*']
    },
    {
        'type': 'nonsense_mutation',
        'rule_name': 'stop_codon_star',
        'regex': f'({aa})(\d+)(\*)',
        'apply_syntax': lambda g: ''.join(g[:2]).upper() + '*',
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'peptide'),
        'format_for_transvar': lambda g, gg: [f'p.{g[0]}{g[1]}*']
    },
    {
        'type': 'partial_amino_acid_deletion',
        'rule_name': 'multiple_aa',
        'regex': f'(?<!{aa})(\d+)\s*[-–]\s*(\d+)(?!{aa})(?:\s+Δaa)?',
        'apply_syntax': lambda g: '-'.join(sorted(g, key=int)).upper(),
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups, gene, 'peptide'),
        'format_for_transvar': lambda g, gg: [f'p.{g[0]}_{g[1]}del'],
    },
    {
        'type': 'partial_amino_acid_deletion',
        'rule_name': 'single_aa',
        'regex': f'(?<!{aa})(\d+)(?!{aa})(?:\s+Δaa)?',
        'apply_syntax': lambda g: g[0],
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups, gene, 'peptide'),
        'format_for_transvar': lambda g, gg: [f'p.{g[0]}del'],
    },
    {
        'type': 'amino_acid_mutation',
        'rule_name': 'CTD',
        'regex': f'(CTD-(?:{ctd_mutation_regex},?\s?)+)$',
        'apply_syntax': lambda g: ctd_apply_syntax(g[0]),
        'further_check': ctd_further_check,
        'format_for_transvar': ctd_format_for_transvar
    },
    {
        'type': 'partial_amino_acid_deletion',
        'rule_name': 'CTD',
        'regex': f'(CTD-(?:{ctd_deletion_regex},?\s?)+)$',
        'apply_syntax': lambda g: ctd_apply_syntax(g[0]),
        'further_check': ctd_further_check,
        'format_for_transvar': ctd_format_for_transvar
    },
    {
        'type': 'amino_acid_deletion_and_mutation',
        'rule_name': 'CTD',
        'regex': f'(CTD-(?:(?:{ctd_mutation_regex}|{ctd_deletion_regex}),?\s?)+)$',
        'apply_syntax': lambda g: ctd_apply_syntax(g[0]),
        'further_check': ctd_further_check,
        'format_for_transvar': ctd_format_for_transvar
    }

]


# Cases contemplated in the below regex:
# AA-23-TT (positive number correctly formatted)
# AA-(-23)-TT (negative number correctly formatted)
# AA23TT (positive number without dashes)
# AA-23TT (negative number without dashes nor parenthesis)
# AA(-23)TT (negative number without dashes)
# AA--23-TT (negative number without parenthesis)
# Note the use of positive and negative lookahead / lookbehind for dashes to include both cases
multi_nt_regex = f'({nt}+)-?((?<=-)(?:-?\d+|\(-\d+\))(?=-)|(?<!-)(?:-?\d+|\(-\d+\))(?!-))-?({nt}+)(?=\\b)'
multi_nt_apply_syntax = lambda g: (''.join(format_negatives(g, [1]))).upper().replace('U', 'T')
multi_nt_check_sequence = lambda g, gg: check_sequence_multiple_pos(g, gg, 'dna')


def multi_nt_format_for_transvar(g: list[str], gene: dict) -> list[str]:
    # There might be parenthesis around negative numbers
    gene_pos = int(g[1].replace('(', '').replace(')', ''))
    # First aminoacid + first position
    first_genome_pos, strand = gene_coords2genome_coords(gene_pos, gene)
    last_genome_pos, _ = gene_coords2genome_coords(gene_pos + len(g[0]) - 1, gene)
    if strand == 1:
        return [f'g.{first_genome_pos}_{last_genome_pos}del{g[0]}ins{g[2]}']
    return [f'g.{last_genome_pos}_{first_genome_pos}del{reverse_complement(g[0])}ins{reverse_complement(g[2])}']


def single_nt_format_for_transvar(g: list[str], gene: dict) -> list[str]:
    # There might be parenthesis around negative numbers
    gene_pos = int(g[1].replace('(', '').replace(')', ''))
    # First aminoacid + first position
    genome_pos, strand = gene_coords2genome_coords(gene_pos, gene)

    if strand == 1:
        return [f'g.{genome_pos}{g[0]}>{g[2]}']
    return [f'g.{genome_pos}{reverse_complement(g[0])}>{reverse_complement(g[2])}']


def insertion_nt_format_for_transvar(g: list[str], gene: dict) -> list[str]:
    # Last nt + last position (it could be written AV23AVLLLL, not ideal maybe, but it passes the regex)
    gene_last_pos = int(g[1].replace('(', '').replace(')', '')) + len(g[0]) - 1
    genome_pos, strand = gene_coords2genome_coords(gene_last_pos, gene)

    if strand == 1:
        return [f'g.{genome_pos}_{genome_pos+1}ins{g[2][1:]}']
    return [f'g.{genome_pos - 1}_{genome_pos}ins{reverse_complement(g[2][1:])}']


def deletion_nt_format_single_for_transvar(g: list[str], gene: dict) -> list[str]:
    gene_pos = int(g[0].replace('(', '').replace(')', ''))
    genome_pos, strand = gene_coords2genome_coords(gene_pos, gene)
    return [f'g.{genome_pos}del']


def deletion_nt_format_multi_for_transvar(g: list[str], gene: dict) -> list[str]:
    gene_pos_start = int(g[0].replace('(', '').replace(')', ''))
    gene_pos_end = int(g[1].replace('(', '').replace(')', ''))
    genome_pos_start, strand = gene_coords2genome_coords(gene_pos_start, gene)
    genome_pos_end, _ = gene_coords2genome_coords(gene_pos_end, gene)

    if strand == 1:
        return [f'g.{genome_pos_start}_{genome_pos_end}del']
    return [f'g.{genome_pos_end}_{genome_pos_start}del']


nucleotide_grammar = [
    {
        'type': 'nucleotide_mutation',
        'rule_name': 'single_nt',
        # Negative numbers are common
        'regex': f'(?<=\\b)({nt}){num}({nt})(?=\\b)',
        'apply_syntax': lambda g: ''.join(format_negatives(g, [1])).upper().replace('U', 'T'),
        'check_sequence': lambda g, gg: check_sequence_single_pos(g, gg, 'dna'),
        'format_for_transvar': single_nt_format_for_transvar,
    },
    {
        'type': 'nucleotide_mutation',
        'rule_name': 'multiple_nt',
        'regex': multi_nt_regex,
        'apply_syntax': multi_nt_apply_syntax,
        'check_sequence': multi_nt_check_sequence,
        # It is only a mutation if the number of nts before and after is the same
        'further_check': lambda g, gg: (len(g[0]) == len(g[2])) & (g[0] != g[2]),
        'format_for_transvar': multi_nt_format_for_transvar,
    },
    {
        'type': 'nucleotide_deletion_and_mutation',
        'rule_name': 'multiple_nt',
        'regex': multi_nt_regex,
        'apply_syntax': multi_nt_apply_syntax,
        'check_sequence': multi_nt_check_sequence,
        # TODO: Here we could even check that a partial deletion has not been written using this syntax: e.g. AVTGLA123AA, but probably rare enough.
        'further_check': lambda g, gg: len(g[0]) > len(g[2]),
        'format_for_transvar': multi_nt_format_for_transvar,
    },
    {
        'type': 'nucleotide_insertion_and_mutation',
        'rule_name': 'multiple_nt',
        'regex': multi_nt_regex,
        'apply_syntax': multi_nt_apply_syntax,
        'check_sequence': multi_nt_check_sequence,
        # We don't want to account insertions here
        'further_check': lambda g, gg: (len(g[0]) < len(g[2])) & (not g[2].startswith(g[0])),
        'format_for_transvar': multi_nt_format_for_transvar,
    },
    {
        'type': 'nucleotide_insertion',
        'rule_name': 'standard',
        'regex': multi_nt_regex,
        'apply_syntax': multi_nt_apply_syntax,
        'check_sequence': multi_nt_check_sequence,
        'further_check': lambda g, gg: (len(g[0]) < len(g[2])) & (g[2].startswith(g[0])),
        'format_for_transvar': insertion_nt_format_for_transvar,
    },
    {
        'type': 'partial_nucleotide_deletion',
        'rule_name': 'usual',
        'regex': f'(?<!{nt}){num}\s*[-–]\s*{num}(?!{nt})',
        'apply_syntax': lambda g: '-'.join(format_negatives(sorted(g, key=lambda x: int(x.replace('(', '').replace(')', ''))), [0, 1])).upper(),
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups, gene, 'dna'),
        'format_for_transvar': deletion_nt_format_multi_for_transvar,
    },
    {
        'type': 'partial_nucleotide_deletion',
        'rule_name': 'single_nt',
        'regex': f'(?<!{nt}){num}(?!{nt})',
        'apply_syntax': lambda g: format_negatives(g, [0])[0],
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups[:1], gene, 'dna'),
        'format_for_transvar': deletion_nt_format_single_for_transvar,
    },
]

# Transition grammars ==================================================

# This grammar recognises the old syntax, and apply_syntax applies the new style
transition_old2new_aminoacid_grammar = copy.deepcopy(aminoacid_grammar_old)

for rule in transition_old2new_aminoacid_grammar:
    if rule['type'] == 'amino_acid_mutation' and rule['rule_name'] == 'multiple_aa':
        rule['apply_syntax'] = lambda g: ''.join(g).upper()
    elif rule['type'] == 'amino_acid_insertion':
        rule['apply_syntax'] = lambda g: f'{g[0]}{g[1]}{g[0]}{g[2]}'.upper()


# Same for nucleotides
transition_old2new_nucleotide_grammar = copy.deepcopy(nucleotide_grammar_old)
for rule in transition_old2new_nucleotide_grammar:
    if rule['type'] == 'nucleotide_mutation' and rule['rule_name'] == 'multiple_nt':
        rule['apply_syntax'] = lambda g: (''.join(format_negatives(g, [1]))).upper().replace('U', 'T')
    elif rule['type'] == 'nucleotide_insertion':
        rule['apply_syntax'] = lambda g: f'{g[0]}{format_negatives(g[1:2],[0])[0]}{g[0]}{g[2]}'.upper().replace('U', 'T')


transition_new2old_aminoacid_grammar = copy.deepcopy(aminoacid_grammar)
for rule in transition_new2old_aminoacid_grammar:
    if rule['rule_name'] == 'multiple_aa':
        rule['apply_syntax'] = lambda g: '-'.join(g).upper()
    elif rule['type'] == 'amino_acid_insertion':
        rule['apply_syntax'] = lambda g: f'{g[0]}{g[1]}-{g[2][1:]}'.upper()

transition_new2old_nucleotide_grammar = copy.deepcopy(nucleotide_grammar)
for rule in transition_new2old_nucleotide_grammar:
    if rule['rule_name'] == 'multiple_nt':
        rule['apply_syntax'] = lambda g: ('-'.join(format_negatives(g, [1]))).upper().replace('U', 'T')
    elif (rule['type'] == 'nucleotide_insertion'):
        rule['apply_syntax'] = lambda g: f'{g[0]}{format_negatives(g[1:2],[0])[0]}-{g[2][1:]}'.upper().replace('U', 'T')

