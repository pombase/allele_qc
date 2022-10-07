from Bio.SeqFeature import FeatureLocation


def get_nt_at_genome_position(pos: int, gene: dict, contig):
    return contig[get_genome_pos(pos, gene)]


def get_genome_pos(pos: int, gene: dict) -> str:
    loc: FeatureLocation
    if 'CDS' in gene:
        loc = gene['CDS'].location
    else:
        if len(gene) != 2:
            # Error, we cannot read this position
            raise ValueError('cannot read sequence, alternative splicing?')
        # The key is the one that is not 'contig'
        key = next(k for k in gene if k != 'contig')
        loc = gene[key].location

    if loc.strand == 1:
        pos = (loc.start - 1) + 1 * loc.strand
    else:
        pos = (loc.end - 1) + 1 * loc.strand

    return pos
