import pandas
import pickle
from genome_functions import handle_systematic_id_for_allele_qc
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def transvar_variant_to_substitution_dict(variant: str, seq_length: int) -> dict:
    """
    A subsitution dict contains a slice with the indexes of residues that will be replaced
    and the residues they will be replaced by
    """

    if variant in ['.', 'p.(=)']:
        return None

    single_sub = re.match(r'^p\.[A-Z](\d+)([A-Z\*])$', variant)
    if single_sub is not None:
        coord, replace_by = single_sub.groups()
        coord = int(coord) - 1
        if replace_by == '*':
            return {'range': slice(coord, seq_length), 'replace_by': '*'}
        return {'range': slice(coord, coord + 1), 'replace_by': replace_by}

    partial_deletion = re.match(r'^p\.[A-Z](\d+)_[A-Z\*](\d+)del(?:\d+|[A-Z\*]+)$', variant)
    if partial_deletion is not None:
        start, end = partial_deletion.groups()
        start = int(start) - 1
        # Both deleted residues included, so we do not subtract 1 from end
        end = int(end)
        return {'range': slice(start, end), 'replace_by': ''}

    single_deletion = re.match(r'^p\.[A-Z\*](\d+)del[A-Z\*]$', variant)
    if single_deletion is not None:
        start, = single_deletion.groups()
        start = int(start) - 1
        end = start + 1
        return {'range': slice(start, end), 'replace_by': ''}

    insertion = re.match(r'^p\.([A-Z])(\d+)_[A-Z\*](\d+)ins([A-Z]+)$', variant)
    if insertion is not None:
        residue, start, end, insert = insertion.groups()
        start = int(start) - 1
        end = int(end) - 1
        if end - start != 1:
            raise ValueError('Length of insertion is not 1 in', variant)
        # Note how this is an empty slice, so the residue at start will not be replaced
        return {'range': slice(end, end), 'replace_by': insert}

    duplication = re.match(r'^p\.([A-Z])(\d+)_[A-Z\*](\d+)dup([A-Z]+)$', variant)
    if duplication is not None:
        residue, start, end, insert = duplication.groups()
        start = int(start) - 1
        end = int(end)
        if end - start != len(insert):
            raise ValueError('Length of duplication does not match duplicated sequence in', variant)
        # Note how this is an empty slice, so the residue at start will not be replaced
        return {'range': slice(start, start), 'replace_by': insert}

    delins = re.match(r'^p\.[A-Z](\d+)_[A-Z\*](\d+)delins([A-Z\*]+)$', variant)
    if delins is not None:
        start, end, insert = delins.groups()
        start = int(start) - 1
        # Here we don't subtract 1 from end because we want to replace the last residue as well
        end = int(end)
        return {'range': slice(start, end), 'replace_by': insert}

    if 'fs' in variant:
        print('Skipping frameshift variant', variant)
        return None

    raise ValueError(f'variant "{variant}" is not supported')


def variant_sequence_from_subsitution_dicts(sequence: str, substitution_dicts: list[dict]) -> str:

    sequence = sequence.upper()
    # First we apply the substitutions that do not change the length of the sequence
    remainning_substitution_dicts = []

    for substitution_dict in substitution_dicts:
        if len(substitution_dict['replace_by']) == (substitution_dict['range'].stop - substitution_dict['range'].start):
            sequence = sequence[:substitution_dict['range'].start] + substitution_dict['replace_by'] + sequence[substitution_dict['range'].stop:]
        else:
            remainning_substitution_dicts.append(substitution_dict)

    # This is a bit of a hacky trick in which we first substitute the residues that are going to be replaced by a character (first numbers, then
    # lowercase letters), for instance, if the substitution dict is {'range': slice(1, 3), 'replace_by': 'CCCCC'}, and we want to apply it on a sequence
    # 'AAAAA', we first replace the residues that will be replaced by a character called replacement_key, in this case 1, `AAAAA` -> `A11AA`, then
    # we substitute the subsequent occurrences of that replacement_key by CCCCC.
    # The problem with this approach is handling cases where there are substitutions and insertions after a given position, that's why we handle it on two
    # steps.

    all_replacement_keys = [str(i) for i in range(10)] + [chr(i) for i in range(ord('a'), ord('z') + 1)]
    all_replacement_keys.reverse()
    sequence_replacement_dict = dict()
    remainning_substitution_dicts2 = list()

    for substitution_dict in remainning_substitution_dicts:
        substitution_target_length = substitution_dict['range'].stop - substitution_dict['range'].start
        if substitution_target_length == 0:
            remainning_substitution_dicts2.append(substitution_dict)
        else:
            try:
                replacement_key = all_replacement_keys.pop()
            except IndexError:
                raise ValueError('More substitutions than replacement keys')
            sequence = sequence[:substitution_dict['range'].start] + replacement_key * substitution_target_length + sequence[substitution_dict['range'].stop:]
            sequence_replacement_dict[replacement_key] = substitution_dict['replace_by']

    # We sort the insertions in reverse order of start, so that we can insert them in the sequence without having to worry about the shifting
    # of the coordinates
    remainning_substitution_dicts2.sort(key=lambda x: x['range'].start, reverse=True)

    for substitution_dict in remainning_substitution_dicts2:
        # We can insert directly here
        sequence = sequence[:substitution_dict['range'].start] + substitution_dict['replace_by'] + sequence[substitution_dict['range'].stop:]

    # Now we substitute the replacement keys by the actual residues
    for replacement_key, replace_by in sequence_replacement_dict.items():
        sequence = re.sub(replacement_key + '+', replace_by, sequence)

    if not sequence.endswith('*'):
        sequence = sequence + '*'  # Add stop codon if it's not there

    return sequence


def process_row(row, genome):

    allele_qc_id = handle_systematic_id_for_allele_qc(row['systematic_id'], row['allele_name'], genome)

    gene = genome[allele_qc_id]
    if 'CDS' not in gene:
        raise ValueError(f'gene {row["systematic_id"]} does not have a CDS, but has transvar protein coordinates {row["transvar_coordinates"]}')

    peptide_length = len(gene['peptide'])
    substitution_dicts = [transvar_variant_to_substitution_dict(v.split('/')[2], peptide_length) for v in row['transvar_coordinates'].split('|')]
    # Remove the None values
    substitution_dicts = [s for s in substitution_dicts if s is not None]
    return variant_sequence_from_subsitution_dicts(str(gene['peptide']), substitution_dicts)


def main():

    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    allele_data = pandas.read_csv('results/allele_results_transvar.tsv', delimiter='\t', na_filter=False)

    # We keep all protein variants (even if they were not described at the protein level)
    allele_data = allele_data[allele_data['transvar_coordinates'].str.contains('/p.')].copy()
    allele_data['variant_sequence'] = allele_data.apply(lambda row: process_row(row, genome), axis=1)
    with open('results/all_protein_variant_sequences.fasta', 'w') as out_file:
        for i, row in allele_data.iterrows():
            sequence_name = f'{row["systematic_id"]}|{row["allele_name"]}'
            sequence_description = f'{row["allele_type"]}|{row["allele_description"]}'
            # We trim the stop codon
            seq_record = SeqRecord(Seq(row['variant_sequence']), id=sequence_name, description=sequence_description)
            SeqIO.write(seq_record, out_file, 'fasta')


if __name__ == '__main__':
    main()
