#%%
import pickle
from genome_functions import get_feature_location_from_string

get_feature_location_from_string('3236617-3237051')

#%%

from Bio import SeqIO

for seq in SeqIO.parse('data/dummy_contig.contig','embl'):
    pass

# %%

with open('data/genome.pickle', 'rb') as ins:
    contig_genome = pickle.load(ins)

contig_genome['']