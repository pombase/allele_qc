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


def expand_CTD_abbreviations(sequence_position: str) -> str:
    """Expand CTD abbreviations to all positions"""

    abbreviations = {
        "CTD_S2": "S1559,S1566,S1579,S1586,S1593,S1600,S1607,S1614,S1621,S1628,S1635,S1642,S1649,S1656,S1663,S1670,S1677,S1684,S1691,S1698,S1705,S1712,S1719,S1726,S1733,S1740,S1747",
        "CTD_T4": "T1554,T1567,T1581,T1588,T1595,T1602,T1609,T1616,T1623,T1630,T1637,T1644,T1651,T1658,T1665,T1672,T1679,T1686,T1693,T1700,T1707,T1714,T1721,T1728,T1735,T1742,T1749",
        "CTD_S5": "S1555,S1562,S1568,S1575,S1582,S1589,S1596,S1603,S1610,S1617,S1624,S1631,S1638,S1645,S1652,S1659,S1666,S1673,S1680,S1687,S1694,S1701,S1708,S1715,S1722,S1729,S1736,S1743,S1750",
        "CTD_S7": "S1557,S1577,S1584,S1591,S1598,S1605,S1612,S1619,S1626,S1633,S1640,S1647,S1654,S1661,S1668,S1675,S1682,S1689,S1696,S1703,S1710,S1717,S1724,S1731,S1738,S1745,S1752"
    }
    for key in abbreviations:
        sequence_position = sequence_position.replace(key, abbreviations[key])
    return sequence_position


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

    # Expand CTD abbreviations
    data['exploded_sequence_position'] = data['sequence_position'].apply(expand_CTD_abbreviations)

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

