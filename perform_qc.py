from models import SyntaxRule
from refinement_functions import check_allele_description
from grammar import allowed_types, aminoacid_grammar, nucleotide_grammar
import pickle
from load_sequences import fasta_genome
import sys
import pandas

# TODO needs fixing
# Build a dictionary PMID - curs
# pmid2curs_dict = dict()
# with open('data/pubs_and_session_ids.csv') as ins:
#     for line in ins:
#         pmid, curs = line.strip().split(',')
#         pmid2curs_dict[pmid] = curs

with open('data/genome.pickle', 'rb') as ins:
    contig_genome = pickle.load(ins)


def get_invalid_CDS_dict():
    """
    Return the dictionary with CDS error
    """
    return {
        'allele_parts': '',
        'needs_fixing': True,
        'rename_to': '',
        'rules_applied': '',
        'pattern_error': '',
        'invalid_error': 'several transcript or CDS missing',
        'sequence_error': '',
        'change_type_to': ''
    }


def get_invalid_dnaseq_dict():
    """
    Return the dictionary with CDS error
    """
    return {
        'allele_parts': '',
        'needs_fixing': True,
        'rename_to': '',
        'rules_applied': '',
        'pattern_error': '',
        'invalid_error': 'sequence missing',
        'sequence_error': '',
        'change_type_to': ''
    }


def main(input_file: str):
    dict_list = list()
    syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
    syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]
    with open(input_file) as ins:
        ins.readline()
        for line in ins:
            systematic_id, allele_description, gene_name, allele_name, allele_synonym, allele_type, pmid = line.strip().split('\t')
            base_dict = {'systematic_id': systematic_id, 'gene_name': gene_name, 'allele_name': allele_name, 'allele_type': allele_type, 'allele_description': allele_description}
            if 'nucleotide' not in allele_type:
                if systematic_id not in fasta_genome or 'peptide' not in fasta_genome[systematic_id]:
                    dict_list.append(base_dict | get_invalid_CDS_dict())
                else:
                    dict_list.append(base_dict | check_allele_description(allele_description, syntax_rules_aminoacids, allele_type, allowed_types, fasta_genome[systematic_id]))
            else:
                if systematic_id not in contig_genome:
                    dict_list.append(base_dict | get_invalid_dnaseq_dict())
                else:
                    dict_list.append(base_dict | check_allele_description(allele_description, syntax_rules_nucleotides, allele_type, allowed_types, contig_genome[systematic_id]))

    data = pandas.DataFrame.from_records(dict_list)
    data.to_csv('results/allele_results.tsv', sep='\t', index=False)
    data[data['needs_fixing'] == True].to_csv('results/allele_errors.tsv', sep='\t', index=False)
    data[data['needs_fixing'] == True][['allele_description', 'rename_to']].to_csv('results/allele_errors_summarised.tsv', sep='\t', index=False)


if __name__ == "__main__":
    main(sys.argv[1])
