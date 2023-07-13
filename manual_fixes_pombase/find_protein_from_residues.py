import pickle
import json

residues = "S103,S107,S110,S12,S132,S134,S143,S145,S148,S149,S15,S151,S154,S157,S18,S333,S335,S345,S347,S552,S554,S590,T106,T144,T152,T339,Y343"
# residues = "S200,S202,S202,S212,S212,S224,S229,S232,S237,S239,S239,S309,S309,S316,S337,S337,S345,S345,S354,S354,S376,S381,S381,S383,S383,S4,S4,S409,S411,S419,S426,S434,S445,S447,S455,S6,T129,T129,T257,T257,T352,T352,T375,T375,T439"
residues = [(r[0], int(r[1:]) - 1) for r in residues.split(",")]
residues = sorted(set(residues), key=lambda x: x[1])
print(residues)

st_the_same = False


def compare_residues(peptide, residue):
    if (residue[0] not in ['S', 'T']) or (not st_the_same):
        return peptide[residue[1]] != residue[0]
    return peptide[residue[1]] not in ['S', 'T']

# Check if the residues are present in any of the sequences (present or past)
with open('../data/genome.pickle', 'rb') as ins:
    genome = pickle.load(ins)

for systematic_id in genome:
    gene = genome[systematic_id]
    if 'peptide' in gene:
        peptide = gene['peptide']
        for residue in residues:
            if (len(peptide) > residue[1]) and compare_residues(peptide, residue):
                break
        else:
            print(systematic_id, 'in current genome')


with open('../data/coordinate_changes_dict.json') as ins:
    coordinate_changes_dict = json.load(ins)

for systematic_id in coordinate_changes_dict:
    for revision in coordinate_changes_dict[systematic_id]:
        peptide = revision['old_alignment'].replace('-', '')
        for residue in residues:
            if (len(peptide) > residue[1]) and compare_residues(peptide, residue):
                break
        else:
            print(systematic_id, 'in past genome', revision['revision'])