import pandas
import pickle
import argparse
from grammar import aminoacid_grammar, nucleotide_grammar
from models import SyntaxRule, find_rule
from transvar_functions import parse_transvar_string, get_transvar_str_annotation, get_anno_db, TransvarAnnotation
from genome_functions import handle_systematic_id_for_qc
from tqdm import tqdm

tqdm.pandas()


def format_for_transvar(row, genome, syntax_rules_aminoacids, syntax_rules_nucleotides) -> str:

    # Transvar uses only gene_ids, while the allele_qc uses a mix to handle multi-transcripts
    gene_systematic_id = row['systematic_id']
    allele_qc_id = handle_systematic_id_for_qc(row, genome)
    gene = genome[allele_qc_id]
    chromosome = gene['contig'].id

    if 'amino_acid' in row['allele_type']:
        syntax_rule = find_rule(syntax_rules_aminoacids, *row['rules_applied'].split(':'))
        prefix = gene_systematic_id
    else:
        syntax_rule = find_rule(syntax_rules_nucleotides, *row['rules_applied'].split(':'))
        prefix = chromosome

    capture_groups = syntax_rule.get_groups(row['allele_parts'])
    transvar_input = '{}:{}'.format(prefix, syntax_rule.format_for_transvar(capture_groups, gene))

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

    # Some variants fall out of the transcript, so we just return the genomic variant
    return '{}/./.'.format(annotations[0].coordinates.split('/')[0])


def get_transvar_coordinates(row, db, genome, exclude_transcripts):
    # print(row['systematic_id'], '<<<>>>', row['transvar_input'])
    allele_qc_id = handle_systematic_id_for_qc(row, genome)
    transcript_id = None if (allele_qc_id == row['systematic_id']) else allele_qc_id
    try:
        transvar_annotation_list = parse_transvar_string(get_transvar_str_annotation('panno' if 'amino_acid' in row['allele_type'] else 'ganno', row['transvar_input'], db))
        return get_transvar_annotation_coordinates(transvar_annotation_list, row['systematic_id'], transcript_id)
    except ValueError as e:
        if e.args[0] == 'no_valid_transcript_found' and row['systematic_id'] in exclude_transcripts:
            return ''
        else:
            raise e


if __name__ == '__main__':
    class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
    parser.add_argument('--genome', default='data/genome.pickle', help='genome dictionary built from contig files.')
    parser.add_argument('--allele_results', default='results/allele_results.tsv')
    parser.add_argument('--exclude_transcripts', default='data/frame_shifted_transcripts.tsv')
    parser.add_argument('--output', default='results/allele_results_transvar.tsv')

    args = parser.parse_args()

    with open(args.genome, 'rb') as ins:
        genome = pickle.load(ins)

    with open(args.exclude_transcripts) as ins:
        exclude_transcripts = set(map(str.strip, ins.readlines()))

    syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
    syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]

    data = pandas.read_csv(args.allele_results, sep='\t', na_filter=False)

    # Remove all errors
    data = data[~data['needs_fixing']].copy()

    # Only allele types that can be fixed
    data = data[data['allele_type'].str.contains('amino_acid|nucleotide', regex=True)].copy()

    # Explode the allele_parts and the rules_applied
    data_exploded = data.copy()
    data_exploded.loc[:, 'allele_parts'] = data_exploded['allele_parts'].apply(str.split, args=['|'])
    data_exploded.loc[:, 'rules_applied'] = data_exploded['rules_applied'].apply(str.split, args=['|'])
    data_exploded = data_exploded.explode(['allele_parts', 'rules_applied'])

    # Apply transvar to each allele_parts
    data_exploded['transvar_input'] = data_exploded.apply(format_for_transvar, axis=1, args=(genome, syntax_rules_aminoacids, syntax_rules_nucleotides))

    anno_db = get_anno_db()
    print('Running transvar on variants... (will take a while)')
    data_exploded['transvar_coordinates'] = data_exploded.progress_apply(get_transvar_coordinates, args=(anno_db, genome, exclude_transcripts), axis=1)

    aggregated_data = data_exploded[['systematic_id', 'allele_description', 'allele_type', 'transvar_coordinates']].groupby(['systematic_id', 'allele_description', 'allele_type'], as_index=False).agg({'transvar_coordinates': '|'.join})

    data.merge(aggregated_data, on=['systematic_id', 'allele_description', 'allele_type'], how='left').to_csv(args.output, sep='\t', index=False)