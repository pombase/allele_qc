
import re


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

allowed_regex = {
    'amino_acid_mutation': [f'{aa}\d+{aa}(,{aa}\d+{aa})*',r'{aa}+-\d+-{aa}+']
}
aa = r'[GPAVLIMCFYWHKRQNEDST]'
def allele_is_valid(allele_type, allele_description):
    match allele_type:
        case 'amino_acid_mutation':
            if re.match('{aa}\d+{aa}(,{aa}\d+{aa})*', allele_description):
                return False
            splitted_syntax_match = re.match('{aa}+-\d+-{aa}', allele_description)
            if splitted_syntax_match:
                print(splitted_syntax_match)
    return False