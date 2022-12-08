import re
from Bio import pairwise2
from genome_functions import get_other_index_from_alignment, get_feature_location_from_string
import pickle
import pandas
from Bio.Seq import Seq
import regex


def old_coords_fix(systematic_id, targets):

    modif_coords = pandas.read_csv('data/only_modified_coordinates.tsv', sep='\t', na_filter=False)

    modif_coords = modif_coords[modif_coords.systematic_id == systematic_id].copy()

    if modif_coords.empty:
        return ''

    # The last coordinate is the last one with 'added'
    last_coord = list(modif_coords['value'])[0]

    prev_coords = list(modif_coords['value'][1:][modif_coords.added_or_removed == 'removed'])
    revisions = list(modif_coords['revision'][1:][modif_coords.added_or_removed == 'removed'])

    with open('data/genome.pickle', 'rb') as ins:
        contig_genome = pickle.load(ins)

    out_str = f'current coordinates: {last_coord}\n'

    for prev_coord, revision in zip(prev_coords, revisions):

        new_feature_loc = get_feature_location_from_string(last_coord)
        old_feature_loc = get_feature_location_from_string(prev_coord)

        new_seq = new_feature_loc.extract(contig_genome[systematic_id]['contig']).translate()
        old_seq = old_feature_loc.extract(contig_genome[systematic_id]['contig']).translate()

        # An alignment where gaps have the maximum penalty, to avoid scattered matches.
        alignments = pairwise2.align.globalms(new_seq.seq, old_seq.seq, match=2, mismatch=-1, open=-100, extend=0, penalize_end_gaps=False)

        new_alignment = alignments[0].seqA
        old_alignment = alignments[0].seqB

        numbers = list()
        values = list()
        for t in targets:
            numbers.append(int(re.search(r'\d+', t).group()) - 1)
            values.append(t[0])
        out_str += '\n'
        out_str += f'coordinates in revision {revision}: {prev_coord}\n'

        for v, old_index in zip(values, numbers):
            new_index = get_other_index_from_alignment(old_alignment, new_alignment, old_index)
            out_str += f'old {v}{old_index+1} > {new_seq[new_index]}{new_index+1}\n'
        out_str += '\n'

    return out_str


def multi_shift_fix(seq, targets):

    numbers = list()
    values = list()
    new_values = list()
    for t in targets:
        numbers.append(int(re.search(r'\d+', t).group()))
        values.append(t[0])
        if not t[-1].isdigit():
            new_values.append(t[-1])
        else:
            new_values.append('')

    numbers_zero = list()
    for n in numbers:
        numbers_zero.append(n - numbers[0])

    i = 0
    out_str = ''
    max_number = len(seq) - max(numbers_zero)
    while i < max_number:
        all_vals = list()
        new_coords = ''
        for n, value, new_value in zip(numbers_zero, values, new_values):
            all_vals.append(value == seq[n + i])
            new_coords += (value + str(n + i + 1) + new_value + ',')
        if all(all_vals):
            out_str += new_coords + '\n'
        i += 1

    return out_str


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
