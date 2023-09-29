from genome_functions import get_nt_at_gene_coord, gene_coords2genome_coords
from Bio.Seq import reverse_complement

# Variable with all the nucleotides, to be used in the regex
# We allow the U, but we will replace it by T
nt = 'ACGUT'
nt = nt + nt.lower()
nt = f'[{nt}]'

# Variable with all aminoacids, to be used in the regex
aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

# We favour formatting negative numbers with parenthesis
# this regex captures both positive and negative numbers without parenthesis, and
# negative numbers with parenthesis
num = '(\(-\d+\)|(?<!\()-?\d+(?!\)))'

# In the multi_aa and multi_nt cases, there are a lot of re-used regex, so we use variables to avoid repetition
multi_aa_regex = f'(?<=\\b)({aa}+)-?(\d+)-?({aa}+)(?=\\b)'
multi_aa_apply_syntax = lambda g: ''.join(g).upper()
multi_aa_check_sequence = lambda g, gg: check_sequence_multiple_pos(g, gg, 'peptide')

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


def format_negatives(input_list: list[str], indexes: list[int]):
    output_list = list(input_list[:])
    for index in indexes:
        this_value = output_list[index].replace('(', '').replace(')', '')
        output_list[index] = this_value if int(this_value) > 0 else f'({this_value})'
    return output_list


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
