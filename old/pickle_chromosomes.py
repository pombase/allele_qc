import pickle
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import argparse
import os

genome: dict[str, dict[str, SeqFeature]] = dict()


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('files', metavar='N', type=str, nargs='+',
                    help='files to be read')
parser.add_argument('--output_folder', default='data/pickled_chromosomes', help='output folder for pickle files')
parser.add_argument('--format', default='embl', help='format of the files to be read (for Biopython)')
args = parser.parse_args()

if not os.path.isdir(args.output_folder):
    os.mkdir(args.output_folder)

for f in args.files:
    file_name = os.path.split(os.path.splitext(f)[0])[-1]
    print(f'pickling {file_name}...')
    iterator = SeqIO.parse(f, args.format)
    contig = next(iterator)
    if next(iterator, None) is not None:
        raise ValueError(f'multiple sequences in file {f}')

    with open(os.path.join(args.output_folder, f'{file_name}.pickle'), 'wb') as out:
        pickle.dump(genome, out, pickle.HIGHEST_PROTOCOL)
