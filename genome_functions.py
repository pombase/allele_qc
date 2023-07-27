from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Seq import reverse_complement
from Bio.GenBank import _FeatureConsumer
from Bio.SeqRecord import SeqRecord


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

    # Passed coordinates are one-based, but only if they are positive
    if pos < 0:
        pos += 1

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


def extract_main_feature_and_strand(gene: dict, downstream: int, upstream: int, get_utrs: bool = True) -> tuple[SeqRecord, int]:
    if 'CDS' not in gene:

        if len(gene) == 2:
            features = list(gene.keys())
            features.remove('contig')
            start_feature = features[0]
            end_feature = features[0]
            strand = gene[start_feature].location.strand
        else:
            raise ValueError('Only supports genes with CDS or RNA genes with a single feature for now')
    else:
        strand = gene["CDS"].location.strand
        start_feature = "5'UTR" if ("5'UTR" in gene and get_utrs) else "CDS"
        end_feature = "3'UTR" if ("3'UTR" in gene and get_utrs) else "CDS"

    if strand == 1:
        end = gene[end_feature].location.end + downstream
        start = gene[start_feature].location.start - upstream
    else:
        end = gene[start_feature].location.end + upstream
        start = gene[end_feature].location.start - downstream

    # Add translation if it exists
    if 'CDS' in gene:
        feat: SeqFeature = gene['CDS']
        feat.qualifiers['translation'] = gene['peptide']

    return gene['contig'][start:end], strand


def process_systematic_id(systematic_id: str, genome: dict, when_several_transcripts: str) -> str:
    """
    If no multiple transcripts exist, return the systematic id, else
    return the transcript id of the longest transcript or the .1 transcript, depending
    on the value of when_several_transcripts
    """

    if systematic_id in genome:
        return systematic_id

    if when_several_transcripts not in ('longest', 'first'):
        return ValueError('when_several_transcripts must be either "longest" or "first"')

    if when_several_transcripts == 'first' and systematic_id + '.1' in genome:
        return systematic_id + '.1'

    # For loci with multiple transcripts, we need to find the one that is the longest
    i = 1
    longest_transcript_systematic_id = None
    longest_transcript_length = 0
    while systematic_id + '.' + str(i) in genome:
        transcript_seq_record, _ = extract_main_feature_and_strand(genome[systematic_id + '.' + str(i)], 0, 0)
        if len(transcript_seq_record) > longest_transcript_length:
            longest_transcript_length = len(transcript_seq_record)
            longest_transcript_systematic_id = systematic_id + '.' + str(i)
        i += 1
    if longest_transcript_systematic_id is None:
        raise ValueError('Systematic id does not exist')
    else:
        return longest_transcript_systematic_id