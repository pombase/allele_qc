import pandas
from allele_fixes import multi_shift_fix
import pickle

with open('data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)

data = pandas.read_csv('results/protein_modification_results_errors_aggregated.tsv', sep='\t', na_filter=False)


def apply_multi_shift_fix(row):
    # We use at least 3 for the multi-shift
    if row['sequence_position'].count(',') < 3:
        return ''
    peptide_seq = genome[row['systematic_id']]['peptide']
    return multi_shift_fix(peptide_seq, row['sequence_position'].split(','))


data['multi_shift_fix'] = data.apply(apply_multi_shift_fix, axis=1)
data.to_csv('results/dummy.tsv', sep='\t', index=False)
