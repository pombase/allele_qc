from refinement_functions import build_regex2modification, allele_is_invalid
from grammar import allowed_types, modifications
import pickle

# Build a dictionary PMID - curs
pmid2curs_dict = dict()
with open('data/pubs_and_session_ids.csv') as ins:
    for line in ins:
        pmid, curs = line.strip().split(',')
        pmid2curs_dict[pmid] = curs

regex2modification = build_regex2modification(modifications)

with open('data/genome.pickle','rb') as ins:
    genome = pickle.load(ins)

with open('data/alleles.tsv') as ins:
    ins.readline()
    for line in ins:
        systematic_id, allele_description, gene_name, allele_name, allele_synonym, allele_type, pmid = line.strip().split('\t')
        if 'nucleotide' in allele_type:
            continue
        if systematic_id not in genome or 'translation' not in genome[systematic_id]:
            reason = 'several transcripts or CDS missing'
        else:
            reason = allele_is_invalid(allele_description, regex2modification, allele_type, allowed_types, genome[systematic_id])

        # curs = '??????'
        # if pmid in pmid2curs_dict:
        #     curs = 'https://curation.pombase.org/pombe/curs/' + pmid2curs_dict[pmid]
        if reason:
            print(allele_name, allele_description, allele_type, reason, sep='\t')