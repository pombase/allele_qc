from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Seq import reverse_complement
from Bio.GenBank import _FeatureConsumer
from Bio.SeqRecord import SeqRecord
import re


def get_nt_at_gene_coord(pos: int, gene: dict, contig):
    # genome_coord is one-based
    genome_coord, strand = gene_coords2genome_coords(pos, gene)
    if strand == 1:
        return contig[genome_coord - 1]
    return reverse_complement(contig[genome_coord - 1])


def get_CDS_or_RNA_feature(gene) -> SeqFeature:
    """
    Gets the main feature of a gene, RNA for a RNA gene, CDS for a protein-coding gene.
    """
    if 'CDS' in gene:
        return gene['CDS']
    else:
        if len(gene) != 2:
            # Error, we cannot read this position
            raise ValueError('Only supports genes with CDS or RNA genes with a single feature for now')
        # The key is the one that is not 'contig'
        key = next(k for k in gene if k != 'contig')
        return gene[key]


def gene_coords2genome_coords(pos: int, gene: dict) -> tuple[int, int]:

    feature = get_CDS_or_RNA_feature(gene)
    loc = feature.location

    # Passed coordinates are one-based, but only if they are positive
    if pos < 0:
        pos += 1

    if loc.strand == 1:
        pos = loc.start + (pos - 1)
    else:
        pos = loc.end - pos

    # We return one-based coordinates as well
    return pos + 1, loc.strand


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


def handle_systematic_id_for_allele_qc(systematic_id, allele_name, genome: dict) -> str:
    """
    Returns the right systematic_id for the allele:
    - For no multi-transcript (row.systematic_id in genome), return row.systematic_id
    - For multi-transcript in which the allele name starts with the primary name + .1, .2, etc, (e.g. zas1.2-V123A) return that transcript (SPBC1198.04c.2).
    - For other multi-transcript genes, return the first transcript, ending in .1 (SPBC1198.04c.1).
    """

    # If it's in the genome dictionary, return it
    if systematic_id in genome:
        return systematic_id

    # Get the first multiple transcript id, if it is a multi-transcript gene
    first_multi_transcript = process_systematic_id(systematic_id, genome, 'first')

    gene_name = get_CDS_or_RNA_feature(genome[first_multi_transcript]).qualifiers['primary_name'][0]

    # If we have reached here, it means that the systematic_id is from a multi-transcript gene
    # If the allele name contains the primary name .1, .2, etc, (e.g. zas1.2) then we pick that transcript (SPBC1198.04c.2).
    # Otherwise, we pick the first transcript

    transcript_id_regex = '^' + gene_name + '\.(\d+)'
    match = re.search(transcript_id_regex, allele_name)
    if match:
        return systematic_id + '.' + match.groups()[0]
    return first_multi_transcript
