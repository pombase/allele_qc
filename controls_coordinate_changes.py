"""
Some controls run on the coordinate change log prior to using it for the analysis
"""
# %%
import pickle
from load_sequences import fasta_genome
import pandas

with open('data/genome.pickle', 'rb') as ins:
    contig_genome = pickle.load(ins)


# %% See if some of the alleles with sequence errors have multiple revisions where their coordinates changed.

# Load alleles with errors
allele_data = pandas.read_csv('results/allele_errors.tsv', delimiter='\t', na_filter=False)
alleles_with_sequence_errors = allele_data[allele_data['sequence_error'] != '']
ids_sequence_errors = set(alleles_with_sequence_errors['systematic_id'])

# Load coordinate changes
coordinate_data = pandas.read_csv('data/all_coordinate_changes_file.tsv', delimiter='\t', na_filter=False)

# We only consider CDS features in genes that have alleles with sequence errors
coordinate_data = coordinate_data[(coordinate_data['feature_type'] == 'CDS') & (coordinate_data['systematic_id'].isin(ids_sequence_errors))]

# See the columns that only differ in value and added_removed > coordinates were modified
d = coordinate_data.drop(columns=['value', 'added_or_removed'])
logi = d.duplicated(keep=False)
coordinate_modifications = coordinate_data[logi]

# Sort the columns to see which ids have been modified more than once
data2sort = coordinate_modifications.copy()
data2sort['count'] = data2sort.groupby('systematic_id')['systematic_id'].transform('count')
sorted_data = data2sort.sort_values(['count', 'systematic_id', 'revision', 'added_or_removed'], ascending=[False, False, False, True])

# Print the data to csv, inspect and see that only nup189 alleles are affected by this problem, and already fixed, by the 1-index fix
sorted_data.to_csv('results/sorted_changes.tsv', sep='\t', index=False)

# %% See if some of the concerned sequences are different in the contigs and the fasta files from pombase
# When I first ran this, the only error was in SPCC162.04c, because it had multiple transcripts

for id in set(sorted_data['systematic_id']):
    if 'peptide' not in contig_genome[id]:
        print(id, 'multipe transcripts')
    elif fasta_genome[id]['peptide'] != contig_genome[id]['peptide']:
        print(id, 'peptides sequences differ')

# %% See if any publication with concerned genes has only deletion alleles (we cannot be sure which sequence reference was used)

all_alleles = pandas.read_csv('data/alleles.tsv', delimiter='\t', na_filter=False)

concerned_alleles = all_alleles[all_alleles['systematic_id'].isin(set(coordinate_data['systematic_id']))]


all_pmids = set()

for pmid_list in concerned_alleles['reference']:
    all_pmids.update(pmid_list.split(','))

pmids_with_mutations = set()
for i, row in concerned_alleles.iterrows():
    if 'mutation' in row['allele_type'] or 'unknown' in row['allele_type']:
        # if row['Allele name'] in alleles_with_sequence_errors['allele_name']:
        pmids_with_mutations.update(row['reference'].split(','))

ambiguous = all_pmids - pmids_with_mutations

ambiguous_data = list()
for i, row in concerned_alleles.iterrows():
    for pmid in row['reference'].split(','):
        if pmid in ambiguous:
            ambiguous_data.append(row)

pandas.DataFrame(ambiguous_data).to_csv('results/ambiguous_coordinate_changes.tsv', sep='\t', index=False)
