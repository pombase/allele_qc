import pandas
from grammar import check_sequence_single_pos, aa
import pickle
import re

with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)

data = pandas.read_csv('data/pombase-chado.modifications', sep='\t', na_filter=False)
# TODO what is the unknown column?
data.columns = ['systematic_id', 'primary_name', 'modification', 'evidence', 'sequence_position', '', 'reference', 'unknown', 'date']

data = data[data['sequence_position'] != '']


def check_func(row):
    if row['systematic_id'] not in genome:
        return 'systematic_id not in genome'

    gene = genome[row['systematic_id']]

    if 'CDS' not in gene:
        return 'gene has no CDS'

    match_obj = re.match(f'({aa})(\d+)', row['sequence_position'])

    if not match_obj:
        return 'pattern error'

    return check_sequence_single_pos(match_obj.groups(), gene, 'peptide')


data['sequence_error'] = data.apply(check_func, axis=1)
data.to_csv('results/protein_modification_results.tsv', sep='\t', index=False)

data[data['sequence_error'] != ''].to_csv('results/protein_modification_results_errors.tsv', sep='\t', index=False)

# Aggregate the errors
