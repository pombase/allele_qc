from Bio.SeqFeature import FeatureLocation
from Bio.Seq import reverse_complement
from Bio.GenBank import _FeatureConsumer


def get_nt_at_genome_position(pos: int, gene: dict, contig):
    genome_coord, strand = gene_coords2genome_coords(pos, gene)
    if strand == 1:
        return contig[genome_coord]
    return reverse_complement(contig[genome_coord])


def gene_coords2genome_coords(pos: int, gene: dict) -> str:
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

    # Passed coordinates are one-based

    if loc.strand == 1:
        pos = loc.start + (pos - 1)
    else:
        pos = loc.end - pos

    return pos, loc.strand


def get_feature_location_from_string(location_str: str) -> FeatureLocation:
    fc = _FeatureConsumer(use_fuzziness=False)
    # We need to initialize a dummy feature
    fc._cur_feature = FeatureLocation(1, 2, 1)
    fc.location(location_str)
    return fc._cur_feature.location


def get_sequence_from_location_string(genome, location_str):
    loc = get_feature_location_from_string(location_str)
    return loc.extract(genome)


def get_other_index_from_alignment(this_alignment, other_alignment, this_index):
    """
    From an alignment of two sequences (this, other), convert the index in
    this to the index in other. Get None if that does not exist. E.g.:

    this_alignment  = "VAQCIKVTVIFLAQCVKVTVIFLAAA"
    other_alignment = "VAQCIKVT----AQCVKVTVIFL"

    this_index == 0  -> returns 0
    this_index == 8  -> returns None (the V does not exist in other_alignment)
    this_index == 12 -> returns 8
    this_index == 19 -> returns None (the A does not exist in other_alignment)

    this_alignment  = "VAQCIKVT----AQCVKVTVIFL"
    other_alignment = "VAQCIKVTVIFLAQCVKVTVIFL"

    this_index == 0  -> returns 0
    this_index == 8  -> returns 12
    """
    count_other = -1
    count_this = -1

    for i in range(len(this_alignment)):
        count_this += this_alignment[i] != '-'
        if i >= len(other_alignment):
            return None
        count_other += other_alignment[i] != '-'
        if count_this == this_index:
            if other_alignment[i] == '-':
                return None
            return count_other
    # The coordinate does not exist in old one (new one is longer)
    return None
