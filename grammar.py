import copy
from ctd_support import ctd_apply_syntax, ctd_further_check, ctd_format_for_transvar
from ctd_support import ctd_mutation_regex, ctd_deletion_regex, ctd_abbreviations
from grammar_funs_and_vars import aa, nt, num, multi_nt_regex, multi_nt_apply_syntax, multi_nt_check_sequence, \
    multi_aa_regex, multi_aa_apply_syntax, multi_aa_check_sequence, \
    multi_aa_format_for_transvar, insertion_aa_format_for_transvar, multi_nt_format_for_transvar, \
    check_sequence_single_pos, check_sequence_multiple_pos, check_multiple_positions_dont_exist, \
    format_negatives, single_nt_format_for_transvar, deletion_nt_format_multi_for_transvar, \
    insertion_nt_format_for_transvar, deletion_nt_format_single_for_transvar


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

disruption_grammar = [
    {
        'type': 'disruption',
        'rule_name': 'usual',
        'regex': '([a-zA-Z]{3}\d+|SP[A-Z0-9]+\.[A-Za-z0-9]+)::(.+?)\+?\s*$',
        'apply_syntax': lambda g: f'{g[0]}::{g[1]}',
    }
]


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
        'regex': f'(?<!{aa}|\d)(\d+)(?!{aa}|\d)(?:\s+Δaa)?',
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
        'regex': f'(CTD-(?:(?:{ctd_mutation_regex}|{ctd_deletion_regex}),?\s?)+|{ctd_abbreviations})$',
        'apply_syntax': lambda g: ctd_apply_syntax(g[0]),
        'further_check': ctd_further_check,
        'format_for_transvar': ctd_format_for_transvar
    }

]

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
        'regex': f'(?<!{nt}|\d){num}(?!{nt}|\d)',
        'apply_syntax': lambda g: format_negatives(g, [0])[0],
        'check_sequence': lambda groups, gene: check_multiple_positions_dont_exist(groups[:1], gene, 'dna'),
        'format_for_transvar': deletion_nt_format_single_for_transvar,
    },
]

## Old grammars (not used anymore) =======

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
        'regex': f'(?<!{aa}|\d)(\d+)(?!{aa}|\d)(?:\s+Δaa)?',
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
        'regex': f'(?<!{nt}|\d){num}(?!{nt}|\d)',
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

