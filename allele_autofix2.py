import pandas
from grammar import aminoacid_grammar
from models import SyntaxRule
from refinement_functions import split_multiple_aa
import pickle
import json
import re
from common_autofix_functions import apply_multi_shift_fix, apply_old_coords_fix, apply_histone_fix, get_preferred_fix, format_auto_fix


with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)

data = pandas.read_csv('results/protein_modification_results_errors_aggregated.tsv', sep='\t', na_filter=False)

with open('results/coordinate_changes_dict2.json') as ins:
    coordinate_changes_dict = json.load(ins)

syntax_rules = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
syntax_rules_dict = {f'{r.type}:{r.rule_name}': r for r in syntax_rules}

data = pandas.read_csv('results/allele_results_errors.tsv', sep='\t', na_filter=False)

# We only want alleles with sequence errors, with mutations described at the aminoacid level
subset_logi = (data['sequence_error'] != '') & (data['rules_applied'].str.contains('amino_acid') | data['rules_applied'].str.contains('nonsense'))
data_subset = data.loc[subset_logi, ['systematic_id', 'allele_description', 'reference', 'change_description_to', 'rules_applied', 'sequence_error']]

# Keep correct description only
data_subset.loc[data_subset['change_description_to'] != '', 'allele_description'] = data_subset.loc[data_subset['change_description_to'] != '', 'change_description_to']

# Explode
cols2explode = ['allele_description', 'rules_applied', 'sequence_error']
for explode_column, sep in zip(cols2explode, [',', '|', '|']):
    data_subset.loc[:, explode_column] = data_subset[explode_column].apply(str.split, args=[sep])
data_subset = data_subset.explode(cols2explode)

# Explode the multiple_aa ones
multi_aa = data_subset['rules_applied'] == 'amino_acid_mutation:multiple_aa'
data_subset.loc[multi_aa, 'allele_description'] = data_subset.loc[multi_aa, 'allele_description'].apply(split_multiple_aa, args=[syntax_rules_dict['amino_acid_mutation:multiple_aa'].regex])
data_subset = data_subset.explode('allele_description')

# Sort the mutations by the first number in them (there might be two numbers in the deletions)
data_subset.loc[:, 'sorting_col'] = data_subset['allele_description'].apply(lambda x: int(re.search(r'\d+', x).group()))
data_subset.sort_values('sorting_col', inplace=True)

# Drop duplicates that may have arised when splitting allele description at either step
data_subset.drop(columns='rules_applied', inplace=True)
data_subset.drop_duplicates(inplace=True)

# Aggregate by publication and systematic id
aggregated_data = data_subset.groupby(['systematic_id', 'reference'], as_index=False).agg({'allele_description': ','.join})


print('applying fixes...')
extra_cols = aggregated_data.apply(apply_old_coords_fix, axis=1, result_type='expand', args=[coordinate_changes_dict, 'allele_description'])
aggregated_data.loc[:, 'old_coords_fix'] = extra_cols.iloc[:, 0]
aggregated_data.loc[:, 'old_coords_revision'] = extra_cols.iloc[:, 1]
aggregated_data.loc[:, 'old_coords_location'] = extra_cols.iloc[:, 2]
aggregated_data['multi_shift_fix'] = aggregated_data.apply(apply_multi_shift_fix, axis=1, args=[genome, 'allele_description'])
aggregated_data['histone_fix'] = aggregated_data.apply(apply_histone_fix, axis=1, args=[genome, 'allele_description'])
aggregated_data.to_csv('a.tsv', sep='\t', index=False)
exit()
extra_cols = aggregated_data.apply(get_preferred_fix, axis=1, result_type='expand')
aggregated_data.loc[:, 'auto_fix_to'] = extra_cols.iloc[:, 0]
aggregated_data.loc[:, 'auto_fix_comment'] = extra_cols.iloc[:, 1]
aggregated_data.rename(columns={'allele_description': 'auto_fix_from'}, inplace=True)
