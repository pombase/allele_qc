from refinement_functions import build_regex2modification, allele_is_invalid
from grammar import allowed_types, modifications


# Build a dictionary PMID - curs
pmid2curs_dict = dict()
with open('data/pubs_and_session_ids.csv') as ins:
    for line in ins:
        pmid, curs = line.strip().split(',')
        pmid2curs_dict[pmid] = curs

regex2modification = build_regex2modification(modifications)

with open('data/alleles.tsv') as ins:
    ins.readline()
    for line in ins:
        systematic_id, allele_description, gene_name, allele_name, allele_synonym, allele_type, pmid = line.strip().split('\t')
        if 'nucleotide' in allele_type:
            continue
        reason = allele_is_invalid(allele_description, regex2modification, allele_type, allowed_types)
        curs = '??????'
        if pmid in pmid2curs_dict:
            curs = 'https://curation.pombase.org/pombe/curs/' + pmid2curs_dict[pmid]
        if reason:
            print(allele_name, allele_description,reason, pmid, curs, sep='\t')