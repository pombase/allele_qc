"""
    amino_acid_deletion_and_mutation
    amino_acid_insertion
    amino_acid_insertion_and_deletion
    amino_acid_insertion_and_mutation
    amino_acid_mutation
    deletion
    disruption
    nonsense_mutation
    nucleotide_insertion
    nucleotide_mutation
    other
    partial_amino_acid_deletion
    partial_nucleotide_deletion
"""
#%%

import re
aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

def allele_is_expected(allele_description, expected):
    if allele_description == expected:
        if allele_description == expected:
            return False
        else:
            return f'typo\tfixto: {expected}'

patterns_allowed = {
    'amino_acid_mutation':
    [
        {
            'regex': f'^({aa})(\d+)({aa})$',
            'apply_syntax': lambda g: ''.join(g).upper(),
            'check_invalid': lambda g: False
        },
        {
            'regex': f'^({aa}+)-?(\d+)-?({aa}+)$',
            'apply_syntax': lambda g: '-'.join(g).upper(),
            'check_invalid': lambda g: f'lengths don\'t match: {g[0]}-{g[2]}' if len(g[0]) != len(g[2]) else False
        },
        {
            'regex': f'^({aa})(\d+)[^a-zA-Z0-9]+(?i:ochre|stop|amber|opal)$',
            'apply_syntax': lambda g: ''.join(g).upper()+'*',
            'check_invalid': lambda g: False
        },
    ]
}

def apply_pattern(match: re.Match, pattern):
    invalid_message = pattern['check_invalid'](match.groups())
    expected = pattern['apply_syntax'](match.groups())
    return invalid_message, expected


def allele_is_invalid(allele_type, allele_description):
    # if ' ' in allele_description:
    #     print(re.match)
    #     return re.sub('\s+','',allele_description)
    # if re.match(r'\s+',allele_description):
    #     return 'allele contains spaces\tfixto: ' + re.sub('\s+','',allele_description)
    allele_parts = allele_description.split(',')
    if allele_type != 'amino_acid_mutation':
        return False
    expected_list = list()
    for part in allele_parts:
        part_matches_regex = False
        for pattern in patterns_allowed['amino_acid_mutation']:
            match = re.match(pattern['regex'], part)
            if match:
                invalid_message, expected = apply_pattern(match, pattern)
                if invalid_message:
                    return invalid_message
                expected_list.append(expected)
                part_matches_regex = True
                break
        if not part_matches_regex:
            return 'pattern not matched: ' + part
    expected_description = ','.join(expected_list)
    if expected_description != allele_description:
        return f'typo\tfixto: {expected_description}'
    return False
    # if re.match(f'{aa}\d+{aa}(,{aa}\d+{aa})*', allele_description):
    #     return False
    # splitted_syntax_match = re.match(f'({aa}+)-?(\d)+-?({aa}+)', allele_description)
    # if splitted_syntax_match:
    #     aa1, number , aa2 = splitted_syntax_match.groups()
    #     if len(aa1) == len(aa2):
    #         return allele_is_expected(allele_description, f'{aa1.upper()}-{number}-{aa2.upper()}')
    #     else:
    #         return 'length of groups differs'
    # elif allele_type == 'amino_acid_insertion':
    #     if re.match(f'{aa}?\d-?{aa}+', allele_type):
    #         return False

    return 'pattern not matched'


with open('data/alleles.tsv') as ins:
    ins.readline()
    for line in ins:
        systematic_id, allele_description, gene_name, allele_name, allele_synonym, allele_type = line.strip().split('\t')
        if allele_type != 'amino_acid_mutation':
            continue
        reason = allele_is_invalid(allele_type, allele_description)
        if reason:
            # print(line.strip(),reason,sep='\t')
            print(allele_description,reason,sep='\t')

