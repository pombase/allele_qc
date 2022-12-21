import pandas
from grammar import check_sequence_single_pos, aa
from refinement_functions import replace_allele_features, sort_result
import pickle

with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)

data = pandas.read_csv('data/pombase-chado.modifications', sep='\t', na_filter=False)
# TODO what is the unknown column?
data.columns = ['systematic_id', 'primary_name', 'modification', 'evidence', 'sequence_position', 'annotation_extension', 'reference', 'unknown', 'date']

data = data[data['sequence_position'] != '']


def check_func(row):
    if row['systematic_id'] not in genome:
        return 'systematic_id not in genome', ''

    gene = genome[row['systematic_id']]

    if 'CDS' not in gene:
        return 'missing_CDS', ''

    result = replace_allele_features([f'(?<!{aa})({aa})(\d+){aa}?'], [row['sequence_position']], [])

    # Extract the matched and unmatched elements
    matches, unmatched = sort_result(result)

    if len(unmatched):
        return 'pattern_error', ''

    correct_name = ','.join(''.join(match.groups()) for match in matches)

    change_sequence_position_to = ''
    if correct_name != row['sequence_position']:
        change_sequence_position_to = correct_name

    errors = [check_sequence_single_pos(match.groups(), gene, 'peptide') for match in matches]
    return '|'.join(errors) if any(errors) else '', change_sequence_position_to


extra_cols = data.apply(check_func, axis=1, result_type='expand')
data.loc[:, 'sequence_error'] = extra_cols.loc[:, 0]
data.loc[:, 'change_sequence_position_to'] = extra_cols.loc[:, 1]
# data.loc[:, ['sequence_error', 'change_sequence_position_to']] = data.apply(check_func, axis=1, result_type='expand')
data.sort_values(['systematic_id', 'sequence_position'], inplace=True)
data.to_csv('results/protein_modification_results.tsv', sep='\t', index=False)

error_data = data[(data['sequence_error'] != '') | (data['change_sequence_position_to'] != '')].copy()
error_data.to_csv('results/protein_modification_results_errors.tsv', sep='\t', index=False)

# Aggregate the errors
sequence_error_data = error_data[~error_data['sequence_error'].isin(['', 'pattern_error', 'missing_CDS'])].copy()
sequence_error_data.loc[sequence_error_data['change_sequence_position_to'] != '', 'sequence_position'] = sequence_error_data['change_sequence_position_to']

sequence_error_data.loc[:, 'sequence_position'] = sequence_error_data['sequence_position'].apply(str.split, args=',')

sequence_error_data = sequence_error_data.explode('sequence_position')
sequence_error_data.loc[:, 'sequence_position'] = sequence_error_data['sequence_position'].astype(str)


sequence_error_data.loc[:, 'sorting_col'] = sequence_error_data['sequence_position'].apply(lambda x: int(x[1:]))
sequence_error_data.sort_values('sorting_col', inplace=True)
aggregated_sequence_error_data = sequence_error_data[['systematic_id', 'reference', 'sequence_position', 'sequence_error']].drop_duplicates().groupby(['systematic_id', 'reference'], as_index=False).agg({'sequence_position': ','.join, 'sequence_error': lambda x: '|'.join(x) if any(x) else ''})
aggregated_sequence_error_data.to_csv('results/protein_modification_results_errors_aggregated.tsv', sep='\t', index=False)
