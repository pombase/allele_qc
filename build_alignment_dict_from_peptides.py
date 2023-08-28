"""
Build a dictionary of alignments based on current and all previous peptide coordinates.

The input is the all_previous_seqs.tsv file generated in https://github.com/pombase/all_previous_sgd_peptide_sequences

The dictionary structure, where keys are the systematic_id of genes:

"SPAC23E2.02": [{
        "new_alignment": "--------------------------------------MNTSENDP ... GYNGTRY*",
        "old_alignment": "MPLGRSSWICCAKYFVNTKSRFNEILPPRFTLIVSFYSMNTSENDP ... SGYNGTRY*"
}],

In the alignment gaps have the maximum penalty, to avoid scattered matches.

-----CAV
MAACATAV

And not

---C--AV
MAACATAV

The reason for this is that the alignment is based on changing / removing introns or changing the start of ending
coordinate of the start or end of the CDS, so you want maximal identity with minimum number of gaps.
"""

import argparse
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json


def main(previous_seqs_file, current_seqs_file, output_file):
    old_seq_dict = dict()
    # This one is a tsv file
    with open(previous_seqs_file) as ins:
        for line in map(str.rstrip, ins):
            ls = line.split('\t')
            gene_id = ls[0]
            if gene_id in old_seq_dict:
                old_seq_dict[ls[0]].append((Seq(ls[1]), ls[2]))
            else:
                old_seq_dict[ls[0]] = [(Seq(ls[1]), ls[2])]

    changes_dict = dict()
    record: SeqRecord
    for record in SeqIO.parse(current_seqs_file, 'fasta'):

        if record.id not in old_seq_dict:
            continue

        changes_dict[record.id] = list()
        for old_seq, revision in old_seq_dict[record.id]:

            alignments = pairwise2.align.globalms(record.seq, old_seq, match=1, mismatch=-2, open=-2, extend=0, penalize_end_gaps=False)
            if len(alignments) == 0:
                print('> No alignment found for {}, skipping'.format(record.id))
                continue
            changes_dict[record.id].append({
                'revision': revision,
                'new_alignment': alignments[0].seqA,
                'old_alignment': alignments[0].seqB
            })

    with open(output_file, 'w') as out:
        json.dump(changes_dict, out, indent=4)


if __name__ == '__main__':
    class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
    parser.add_argument('--previous_seqs', default='data/sgd/all_previous_seqs.tsv')
    parser.add_argument('--current_seqs', default='data/sgd/current_protein_seqs.fasta')
    parser.add_argument('--output', default='data/sgd/coordinate_changes_dict.json')
    args = parser.parse_args()

    main(args.previous_seqs, args.current_seqs, args.output)
