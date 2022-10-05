# %%

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import glob
# The genome contains a dictionary in which the key is the systematic ID,
# and the value is a dictionary of features with the feature_type as key.
import pickle
genome: dict[str, dict[str, SeqFeature]] = dict()

contig_files = glob.glob('data/*.contig')
out = open('data/genome_qc.tsv', 'w')
for f in contig_files:
    out.write('reading: ' + f + '\n')
    iterator = SeqIO.parse(f, 'embl')
    contig = next(iterator)
    if next(iterator, None) is not None:
        raise ValueError(f'multiple sequences in file {f}')
    for feature in contig.features:
        feature: SeqFeature
        if 'systematic_id' not in feature.qualifiers:
            continue
        gene_id = feature.qualifiers['systematic_id'][0]
        if gene_id not in genome:
            genome[gene_id] = dict()
        feature_type = feature.type
        if feature_type in ['intron', 'misc_feature']:
            continue
        if feature_type in genome[gene_id]:
            raise ValueError(
                f'several features of {feature_type} for {gene_id}')

        genome[gene_id][feature_type] = feature
        genome[gene_id]['contig'] = contig
        # if feature_type == 'CDS' and not any([('pseudogene' in prod or 'dubious' in prod) for prod in feature.qualifiers['product']]):
        if feature_type == 'CDS':
            cds_seq = feature.extract(contig).seq
            genome[gene_id]['translation'] = cds_seq.translate()
            errors = list()
            if len(cds_seq) % 3 != 0:
                errors.append('CDS length not multiple of 3')
            if genome[gene_id]['translation'][-1] != '*':
                errors.append('does not end with STOP codon')
            if genome[gene_id]['translation'].count('*') > 1:
                errors.append('multiple stop codons')
            if len(errors):
                out.write(gene_id + '\t' + ','.join(errors) + '\t' + str(feature.qualifiers['product']) + '\n')

        genome['contig'] = contig
out.close()

with open('data/genome.pickle', 'wb') as out:
    pickle.dump(genome, out, pickle.HIGHEST_PROTOCOL)
