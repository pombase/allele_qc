import pickle

with open('data/genome.pickle', 'rb') as ins:
    contig_genome = pickle.load(ins)

with open('data/fasta_genome.pickle', 'rb') as ins:
    fasta_genome = pickle.load(ins)


for id in contig_genome:
    if id not in fasta_genome:
        continue
    if 'peptide' in contig_genome[id]:
        if 'peptide' not in fasta_genome[id]:
            print(id, 'no peptide')
            continue
        if contig_genome[id]['peptide'] != fasta_genome[id]['peptide']:
            print(id)
            print('contig')
            print(contig_genome[id]['peptide'])
            print(''.join(' ' if contig_genome[id]['peptide'][i] == fasta_genome[id]['peptide'][i] else '|' for i in range(len(contig_genome[id]['peptide']))))
            print(fasta_genome[id]['peptide'])
            print('fasta')
            print()
