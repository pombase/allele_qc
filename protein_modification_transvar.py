"""
Uses transvar to represent the modification positions in standard genomic coordinates.

Removes all lines with sequence errors.

"""

import pandas
import pickle
import argparse
from transvar_functions import parse_transvar_string, get_transvar_str_annotation, get_anno_db, TransvarAnnotation
from genome_functions import process_systematic_id
from tqdm import tqdm

tqdm.pandas()


def format_for_transvar(row, genome):

    # Transvar uses only gene_ids, while the pipeline uses a mix to handle multi-transcripts
    gene_systematic_id = row['systematic_id']
    transvar_input = '{}:p.{}'.format(gene_systematic_id, row['exploded_sequence_position'])

    return transvar_input


def get_transvar_annotation_coordinates(annotations: list[TransvarAnnotation], gene_id: str, transcript_id: str) -> TransvarAnnotation:
    # There may be multiple annotations for the same systematic_id if there are multiple
    # transcripts
    for annotation in annotations:
        if annotation.gene == gene_id:
            if transcript_id is None:
                return annotation.coordinates
            elif annotation.transcript.split()[0] == transcript_id:
                return annotation.coordinates

    raise ValueError('Cannot find annotation for {} {}'.format(gene_id, transcript_id))


def get_transvar_coordinates(row, db, genome, exclude_transcripts):

    # print(row['systematic_id'], '<<<>>>', row['transvar_input'])
    qc_id = process_systematic_id(row['systematic_id'], genome, 'first')
    transcript_id = None if (qc_id == row['systematic_id']) else qc_id
    try:
        transvar_annotation_list = parse_transvar_string(get_transvar_str_annotation('panno', row['transvar_input'], db))
        return get_transvar_annotation_coordinates(transvar_annotation_list, row['systematic_id'], transcript_id)
    except ValueError as e:
        if e.args[0] == 'no_valid_transcript_found' and row['systematic_id'] in exclude_transcripts:
            return ''
        else:
            raise e


def main(genome_file, protein_modification_results_file, exclude_transcripts_file, output_file):

    with open(genome_file, 'rb') as ins:
        genome = pickle.load(ins)

    with open(exclude_transcripts_file) as ins:
        exclude_transcripts = set(map(str.strip, ins.readlines()))

    data = pandas.read_csv(protein_modification_results_file, sep='\t', na_filter=False)

    # Remove sequence errors
    data = data[data['sequence_error'] == ''].copy()

    # Use corrected sequence position if available
    data['exploded_sequence_position'] = data['sequence_position']
    data.loc[data['change_sequence_position_to'] != '', 'exploded_sequence_position'] = data['change_sequence_position_to']

    # Explode the sequence_position and the rules_applied
    data_exploded = data[['systematic_id', 'sequence_position', 'exploded_sequence_position']].copy()
    data_exploded.drop_duplicates(inplace=True)
    data_exploded.loc[:, 'exploded_sequence_position'] = data_exploded['exploded_sequence_position'].apply(str.split, args=[','])
    data_exploded = data_exploded.explode(['exploded_sequence_position'])

    # Apply transvar syntax to each sequence position
    data_exploded['transvar_input'] = data_exploded.apply(format_for_transvar, axis=1, args=(genome,))

    anno_db = get_anno_db('data/pombe_genome.gtf.transvardb', 'data/pombe_genome.fa')
    print('Running transvar on protein modifications... (will take a while)')
    data_exploded['transvar_coordinates'] = data_exploded.progress_apply(get_transvar_coordinates, args=(anno_db, genome, exclude_transcripts), axis=1)

    aggregated_data = data_exploded.groupby(['systematic_id', 'sequence_position'], as_index=False).agg({'transvar_coordinates': '|'.join})

    data = data.merge(aggregated_data, on=['systematic_id', 'sequence_position'], how='left')
    data.drop(columns=['exploded_sequence_position'], inplace=True)
    data.to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
    parser.add_argument('--genome', default='data/genome.pickle', help='genome dictionary built from contig files (see load_genome.py).')
    parser.add_argument('--protein_modification_results', default='results/protein_modification_results.tsv', help='output of protein_modification_qc.py')
    parser.add_argument('--exclude_transcripts', default='data/frame_shifted_transcripts.tsv', help='transcripts to exclude from transvar because they are known to be problematic')
    parser.add_argument('--output', default='results/protein_modification_results_transvar.tsv', help='output file')

    args = parser.parse_args()
    main(args.genome, args.protein_modification_results, args.exclude_transcripts, args.output)

