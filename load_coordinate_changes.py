"""
Produces an output file results/alleles_coordinate_change.tsv, with all alleles which has a column named 'uncertain_coordinate_change'.
This column is True if the allele should not be used when fixing from old to new coordinates, because it is uncertain to
which version of the genome the coordinate refers to.

Runs some additional controls as well. It is quite pombase-specific so it uses the default paths, can be adapted if needed.
"""
import pickle
import pandas

with open('data/genome.pickle', 'rb') as ins:
    contig_genome = pickle.load(ins)

with open('data/fasta_genome.pickle', 'rb') as ins:
    fasta_genome = pickle.load(ins)


# %% See if some of the alleles with sequence errors have multiple revisions where their coordinates changed.

# Load alleles with errors
allele_data = pandas.read_csv('results/allele_results_errors.tsv', delimiter='\t', na_filter=False)
alleles_with_sequence_errors = allele_data[allele_data['sequence_error'] != '']
ids_sequence_errors = set(alleles_with_sequence_errors['systematic_id'])

# Load coordinate changes
coordinate_modifications = pandas.read_csv('data/only_modified_coordinates.tsv', delimiter='\t', na_filter=False)

# We only consider CDS features in genes that have alleles with sequence errors
coordinate_modifications = coordinate_modifications[(coordinate_modifications['feature_type'] == 'CDS') & (coordinate_modifications['systematic_id'].isin(ids_sequence_errors))].copy()

# Sort the columns to see which ids have been modified more than once
data2sort = coordinate_modifications.copy()
data2sort['count'] = data2sort.groupby('systematic_id')['systematic_id'].transform('count')
sorted_data = data2sort.sort_values(['count', 'systematic_id', 'revision', 'added_or_removed'], ascending=[False, False, False, True])

# Save the systematic ids that should be excluded because they have been changed more than once
with open('results/systematic_ids_excluded_coordinate_changes.txt', 'w') as out:
    out.write('\n'.join(set(sorted_data['systematic_id'][sorted_data['count'] > 2])))

# %% See if some of the concerned sequences are different in the contigs and the fasta files from pombase
# When I first ran this, the only error was in SPCC162.04c, because it had multiple transcripts

for id in set(sorted_data['systematic_id']):
    if 'peptide' not in contig_genome[id]:
        print(id, 'multipe transcripts')
    elif fasta_genome[id]['peptide'] != contig_genome[id]['peptide']:
        print(id, 'peptides sequences differ')

# %% See if any publication with concerned genes has only deletion alleles (we cannot be sure which sequence reference was used)

all_alleles = pandas.read_csv('data/alleles.tsv', delimiter='\t', na_filter=False)

concerned_alleles = all_alleles[~all_alleles['allele_type'].str.contains('nucl') & all_alleles['systematic_id'].isin(set(coordinate_modifications['systematic_id']))].copy()


all_pmids = set()

for pmid_list in concerned_alleles['reference']:
    all_pmids.update(pmid_list.split(','))

pmids_with_mutations = set()
for i, row in concerned_alleles.iterrows():
    if 'mutation' in row['allele_type'] or 'unknown' in row['allele_type']:
        # if row['Allele name'] in alleles_with_sequence_errors['allele_name']:
        pmids_with_mutations.update(row['reference'].split(','))

ambiguous = all_pmids - pmids_with_mutations

alleles_to_be_fixed = list()
concerned_alleles['uncertain_coordinate_change'] = False
for i, row in concerned_alleles.iterrows():
    for pmid in row['reference'].split(','):
        if pmid in ambiguous:
            concerned_alleles.at[i, 'uncertain_coordinate_change'] = True
            break


concerned_alleles[concerned_alleles['uncertain_coordinate_change'] == False].to_csv('results/alleles_coordinate_change.tsv', sep='\t', index=False)
