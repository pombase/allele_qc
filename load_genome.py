# %%

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, CompoundLocation
import glob
# The genome contains a dictionary in which the key is the systematic ID,
# and the value is a dictionary of features with the feature_type as key.

genome : dict[str, dict[str,SeqFeature]] = dict()

contig_files = glob.glob('data/*.contig')

for f in contig_files[0:1]:
    print('reading',f)
    iterator = SeqIO.parse(f,'embl')
    contig = next(iterator)
    if next(iterator,None) is not None:
        raise ValueError(f'multiple sequences in file {f}')
    for feature in contig.features:
        feature: SeqFeature
        if not 'systematic_id' in feature.qualifiers:
            continue
        gene_id = feature.qualifiers['systematic_id'][0]
        if gene_id not in genome:
            genome[gene_id] = dict()
        feature_type = feature.type
        if feature_type in ['intron','misc_feature']:
            continue
        if feature_type in genome[gene_id]:
            raise ValueError(f'several features of {feature_type} for {gene_id}')

        genome[gene_id][feature_type] = feature
        if feature_type == 'CDS' and not any([('pseudogene' in prod or 'dubious' in prod) for prod in feature.qualifiers['product']]):
            cds_seq = feature.extract(contig).seq
            genome[gene_id]['translation'] = cds_seq.translate()
            if len(cds_seq) % 3 != 0:
                print(gene_id,'CDS length not multiple of 3')
            elif genome[gene_id]['translation'][-1] != '*':
                print(cds_seq)
                print(genome[gene_id]['translation'])
                print(gene_id,'does not end with STOP codon')
        genome['contig'] = contig

# %%




