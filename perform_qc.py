"""
Runs the main analysis pipeline (see readme). The output is stored in the file indicated
in --output, in addition, two extra output files are created, for the default output value:

Only the subset of alleles that needs fixing.
results/allele_results_errors.tsv

Only the subset of alleles that needs fixing, only the columns 'allele_description' and 'rename_to'.
results/allele_results_errors_summarised.tsv
"""

from models import SyntaxRule
from refinement_functions import check_allele_description
from grammar import allowed_types, aminoacid_grammar, nucleotide_grammar
import pickle
import pandas
import argparse


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--genome', default='data/genome.pickle', help='genome dictionary built from contig files.')
parser.add_argument('--fasta_genome', default='data/fasta_genome.pickle', help='genome dictionary built from fasta files.')
parser.add_argument('--alleles', default='data/alleles.tsv')
parser.add_argument('--output', default='results/allele_results.tsv')
args = parser.parse_args()


with open(args.genome, 'rb') as ins:
    contig_genome = pickle.load(ins)

with open(args.fasta_genome, 'rb') as ins:
    fasta_genome = pickle.load(ins)


def invalid_dict():
    """
    Return the dictionary with error
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


allele_data = pandas.read_csv(args.alleles, delimiter='\t', na_filter=False)
syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]

output_data_list = list()
for i, row in allele_data.iterrows():
    if 'nucleotide' not in row['allele_type']:
        if row.systematic_id not in fasta_genome or 'peptide' not in fasta_genome[row.systematic_id]:
            output_data_list.append(dict(row) | invalid_dict())
        else:
            output_data_list.append(dict(row) | check_allele_description(row.allele_description, syntax_rules_aminoacids, row.allele_type, allowed_types, fasta_genome[row.systematic_id]))
    else:
        if row.systematic_id not in contig_genome:
            output_data_list.append(dict(row) | invalid_dict())
        else:
            output_data_list.append(dict(row) | check_allele_description(row.allele_description, syntax_rules_nucleotides, row.allele_type, allowed_types, contig_genome[row.systematic_id]))

output_data = pandas.DataFrame.from_records(output_data_list)
output_data.to_csv(args.output, sep='\t', index=False)

root_output_name = args.output.split('.')[0]

output_data[output_data['needs_fixing'] == True].to_csv(f'{root_output_name}_errors.tsv', sep='\t', index=False)
output_data[output_data['needs_fixing'] == True][['allele_description', 'rename_to']].to_csv(f'{root_output_name}_errors_summarised.tsv', sep='\t', index=False)
