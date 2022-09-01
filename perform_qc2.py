from refinement_functions import build_regex2modification, allele_is_invalid

aa = 'GPAVLIMCFYWHKRQNEDST'
aa = aa + aa.lower()
aa = f'[{aa}]'

allowed_types = {
    frozenset({'amino_acid_mutation'}): 'amino_acid_mutation',
    frozenset({'partial_amino_acid_deletion'}): 'partial_amino_acid_deletion'
}

modifications = [
    {
        'type': 'amino_acid_mutation',
        'regex': f'(?<!{aa})({aa})(\d+)({aa})(?!{aa})',
        'apply_syntax': lambda g: ''.join(g).upper(),
        'check_invalid': lambda g: False
    },
    {
        'type': 'amino_acid_mutation',
        'regex': f'({aa}+)-?(\d+)-?({aa}+)(?!\d)',
        # We fix the case in which dashes are used for a single aa substitution: K-90-R
        'apply_syntax': lambda g: '-'.join(g).upper() if len(g[0])!=1 else ''.join(g).upper(),
        'check_invalid': lambda g: f'lengths don\'t match: {g[0]}-{g[2]}' if len(g[0]) != len(g[2]) else False
    },
    {
        'type': 'amino_acid_mutation',
        'regex': f'({aa})(\d+)[^a-zA-Z0-9]+(?i:ochre|stop|amber|opal)',
        'apply_syntax': lambda g: ''.join(g).upper()+'*',
        'check_invalid': lambda g: False
    },
    {
        'type': 'amino_acid_mutation',
        'regex': f'({aa})(\d+)(\*|stop)',
        'apply_syntax': lambda g: ''.join(g[:2]).upper()+'*',
        'check_invalid': lambda g: False
    },
    {
        'type': 'partial_amino_acid_deletion',
        'regex': f'(?<!{aa})(\d+)-(\d+)(?!{aa})',
        'apply_syntax': lambda g: '-'.join(g).upper(),
        'check_invalid': lambda g: False
    },
]

regex2modification = build_regex2modification(modifications)

with open('data/alleles.tsv') as ins:
    ins.readline()
    for line in ins:
        systematic_id, allele_description, gene_name, allele_name, allele_synonym, allele_type = line.strip().split('\t')
        if allele_type != 'amino_acid_mutation':
            continue
        reason = allele_is_invalid(allele_description, regex2modification, allele_type, allowed_types)
        if reason:
            print(allele_description,reason,sep='\t')