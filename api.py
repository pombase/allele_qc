from fastapi import FastAPI, HTTPException, Query
import json
from starlette.responses import RedirectResponse, PlainTextResponse, FileResponse
from pydantic import BaseModel
import pickle
from grammar import allowed_types, aminoacid_grammar, nucleotide_grammar, disruption_grammar
from models import SyntaxRule, find_rule
from refinement_functions import check_allele_description, split_multiple_aa
from enum import Enum
from allele_fixes import multi_shift_fix, old_coords_fix, primer_mutagenesis as primer_mutagenesis_func
from common_autofix_functions import apply_histone_fix
from protein_modification_qc import check_func as check_modification_description
import re
from typing import Optional
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio import SeqIO
import tempfile
import os
from starlette.background import BackgroundTask


syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]
syntax_rules_disruption = [SyntaxRule.parse_obj(r) for r in disruption_grammar]
multi_aa_rule = find_rule(syntax_rules_aminoacids, 'amino_acid_mutation', 'multiple_aa')


class DNAorProtein(str, Enum):
    protein = 'protein'
    dna = 'dna'


class AlleleType(str, Enum):
    amino_acid_deletion_and_mutation = 'amino_acid_deletion_and_mutation'
    amino_acid_insertion = 'amino_acid_insertion'
    amino_acid_insertion_and_deletion = 'amino_acid_insertion_and_deletion'
    amino_acid_insertion_and_mutation = 'amino_acid_insertion_and_mutation'
    amino_acid_mutation = 'amino_acid_mutation'
    disruption = 'disruption'
    nonsense_mutation = 'nonsense_mutation'
    nucleotide_insertion = 'nucleotide_insertion'
    nucleotide_mutation = 'nucleotide_mutation'
    other = 'other'
    partial_amino_acid_deletion = 'partial_amino_acid_deletion'
    partial_nucleotide_deletion = 'partial_nucleotide_deletion'
    nucleotide_deletion_and_mutation = 'nucleotide_deletion_and_mutation'
    nucleotide_insertion_and_mutation = 'nucleotide_insertion_and_mutation'


class CheckAlleleRequest(BaseModel):
    systematic_id: str
    allele_description: str
    allele_type: AlleleType

    class Config:
        schema_extra = {
            "example": {
                "systematic_id": "SPBC359.03c",
                "allele_description": "V123A,PLR-140-AAA,150-600",
                "allele_type": "partial_amino_acid_deletion"}
        }


class CheckModificationRequest(BaseModel):
    systematic_id: str
    sequence_position: str

    class Config:
        schema_extra = {
            "example": {
                "systematic_id": "SPBC359.03c",
                "sequence_position": "V123; V124,V125",
            }
        }


class MultiShiftRequest(BaseModel):
    systematic_id: str
    targets: str

    class Config:
        schema_extra = {
            "example": {
                "systematic_id": "SPAPB1A10.09",
                "targets": "S123,A124,N125"
            }
        }


class OldCoordsFixRequest(BaseModel):
    systematic_id: str
    targets: str

    class Config:
        schema_extra = {
            "example": {
                "systematic_id": "SPBC1706.01",
                "targets": 'P170A,V223A,F225A,AEY-171-LLL'
            }
        }


class HistoneFixRequest(BaseModel):
    systematic_id: str
    targets: str

    class Config:
        schema_extra = {
            "example": {
                "systematic_id": "SPAC1834.04",
                "targets": 'ART-1-LLL,K9A,K14R,K14A'
            }
        }


class PrimerRequest(BaseModel):
    systematic_id: str
    primer: str
    max_mismatch: int

    class Config:
        schema_extra = {
            "example": {
                "systematic_id": "SPAPB1A10.09",
                "primer": 'TTAGAGGTTATTAATTCCTAAGAAGAAGAAATTTTGG',
                'max_mismatch': 3
            }
        }


class CheckAlleleDescriptionResponse(BaseModel):

    allele_parts: str
    needs_fixing: bool
    change_description_to: str
    rules_applied: str
    pattern_error: str
    invalid_error: str
    sequence_error: str
    change_type_to: str
    user_friendly_fields: Optional['CheckAlleleDescriptionResponse'] = None


class CheckModificationResponse(BaseModel):

    needs_fixing: bool
    sequence_error: str
    change_sequence_position_to: str
    user_friendly_fields: Optional['CheckModificationResponse'] = None


class AlleleFix(BaseModel):

    values: str


class OldCoordsFix(AlleleFix):

    revision: str
    location: str


def get_modification_user_friendly_fields(obj: CheckModificationResponse) -> CheckModificationResponse:
    new_obj = obj.copy()

    if obj.sequence_error:
        positions_dont_exist = list()
        residues_dont_match = list()
        for e in obj.sequence_error.split('|'):
            if e == '':
                continue
            if e[0].isalpha():
                residues_dont_match.append(e)
            else:
                positions_dont_exist.append(e)
        out_str = ''
        if len(positions_dont_exist):
            out_str = out_str + f'The following sequence positions don\'t exist: {",".join(positions_dont_exist)}'
        if len(residues_dont_match):
            if len(positions_dont_exist):
                out_str = out_str + '\n'
            out_str = out_str + f'The following sequence positions don\'t contain the indicated residues: {",".join(residues_dont_match)}'
        new_obj.sequence_error = out_str

    if obj.change_sequence_position_to:
        new_obj.change_sequence_position_to = f'There is a syntax error in the modification description, it should be changed to: {obj.change_sequence_position_to}'

    return new_obj


def get_allele_user_friendly_fields(obj: CheckAlleleDescriptionResponse) -> CheckAlleleDescriptionResponse:
    new_obj = obj.copy()
    if obj.allele_parts:
        new_obj.allele_parts = f'We identified {len(obj.allele_parts.split("|"))} parts in the allele description, which are: {obj.allele_parts.replace("|",", ")}'
    if obj.change_description_to:
        new_obj.change_description_to = f'There is a syntax error in the allele description, it should be changed to: {obj.change_description_to}'
    if obj.rules_applied:
        reformatted_rules = ', '.join(r.split(':')[0] for r in obj.rules_applied.split('|'))
        new_obj.rules_applied = f'We identified that the allele parts correspond to the following types of mutations: {reformatted_rules}'
    if obj.pattern_error:
        new_obj.pattern_error = f'The following parts of the allele description do not follow the existing syntax: {obj.pattern_error}'
    if obj.invalid_error:
        new_obj.invalid_error = f'Invalid error: {obj.invalid_error}'
    if obj.sequence_error:
        positions_dont_exist = list()
        residues_dont_match = list()
        for e in obj.sequence_error.split('|'):
            if e == '':
                continue

            if e[0].isalpha():
                residues_dont_match.append(e)
            else:
                positions_dont_exist.append(e)
        out_str = ''
        if len(positions_dont_exist):
            out_str = out_str + f'The following sequence positions don\'t exist: {",".join(positions_dont_exist)}'
        if len(residues_dont_match):
            if len(positions_dont_exist):
                out_str = out_str + '\n'
            out_str = out_str + f'The following sequence positions don\'t contain the indicated residues: {",".join(residues_dont_match)}'
        new_obj.sequence_error = out_str
    if obj.change_type_to:
        new_obj.change_type_to = f'The indicated allele_type is not correct based on the existing mutations, it should be changed to {obj.change_type_to}'

    return new_obj

def extract_main_feature_and_strand(gene: dict, downstream: int, upstream: int) -> tuple[SeqRecord, int]:
    if 'CDS' not in gene:

        if len(gene) == 2:
            features = list(gene.keys())
            features.remove('contig')
            start_feature = features[0]
            end_feature = features[0]
            strand = gene[start_feature].location.strand
        else:
            raise HTTPException(400, 'Only supports genes with CDS or RNA genes with a single feature for now')
    else:
        strand = gene["CDS"].location.strand
        start_feature = "5'UTR" if "5'UTR" in gene else "CDS"
        end_feature = "3'UTR" if "3'UTR" in gene else "CDS"

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


app = FastAPI()


@ app.get("/")
async def root():
    return RedirectResponse("/docs")


@ app.get("/check_allele", response_model=CheckAlleleDescriptionResponse)
async def check_allele_get(systematic_id: str = Query(example="SPBC359.03c"), allele_description: str = Query(example="V123A,PLR-140-AAA,150-600"), allele_type: AlleleType = Query(example="partial_amino_acid_deletion")):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    if systematic_id not in genome:
        raise HTTPException(404, 'Systematic ID not found in the genome')
    if 'amino' in allele_type:
        response_data = CheckAlleleDescriptionResponse.parse_obj(
            check_allele_description(allele_description, syntax_rules_aminoacids, allele_type, allowed_types, genome[systematic_id])
        )
    else:
        response_data = CheckAlleleDescriptionResponse.parse_obj(
            check_allele_description(allele_description, syntax_rules_nucleotides, allele_type, allowed_types, genome[systematic_id])
        )

    response_data.user_friendly_fields = get_allele_user_friendly_fields(response_data)
    return response_data


@ app.get("/check_modification", response_model=CheckModificationResponse)
async def check_modification(systematic_id: str = Query(example="SPBC359.03c"), sequence_position: str = Query(example="V123; V124,V125")):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    if systematic_id not in genome:
        raise HTTPException(404, 'Systematic ID not found in the genome')
    errors, change_sequence_position_to = check_modification_description({'systematic_id': systematic_id, 'sequence_position': sequence_position}, genome)
    needs_fixing = errors != '' or change_sequence_position_to != ''
    response_data = CheckModificationResponse(sequence_error=errors, change_sequence_position_to=change_sequence_position_to, needs_fixing=needs_fixing)
    response_data.user_friendly_fields = get_modification_user_friendly_fields(response_data)
    return response_data


@ app.post("/primer")
async def primer_mutagenesis(request: PrimerRequest):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    gene = genome[request.systematic_id]
    if 'CDS' in gene:
        has_peptide = True
        seq = gene['CDS'].extract(gene['contig'])
    else:
        has_peptide = False
        if len(gene) != 2:
            # Error, we cannot read this position
            raise ValueError('cannot read sequence, alternative splicing?')
        # The key is the one that is not 'contig'
        key = next(k for k in gene if k != 'contig')
        seq = gene[key].extract(gene['contig'])

    return PlainTextResponse(primer_mutagenesis_func(seq.seq, request.primer, request.max_mismatch, has_peptide))


@ app.post("/multi_shift_fix", response_model=list[AlleleFix])
async def fix_with_multi_shift(request: MultiShiftRequest):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    gene = genome[request.systematic_id]
    targets_1 = request.targets.split(',')
    targets = list()
    for t in targets_1:
        if not re.match(multi_aa_rule.regex, t):
            targets.append(t)
        else:
            targets.extend(split_multiple_aa(t, multi_aa_rule.regex))

    if len([t for t in targets if t[0].isalpha()]) >= 3:
        result = multi_shift_fix(gene['peptide'], targets)
        return [AlleleFix.parse_obj({'values': i}) for i in result]
    else:
        return []


@ app.post("/old_coords_fix", response_model=list[OldCoordsFix])
async def fix_with_old_coords(request: OldCoordsFixRequest):
    with open('data/coordinate_changes_dict.json') as ins:
        coordinate_changes_dict = json.load(ins)
    targets_1 = request.targets.split(',')
    targets = list()
    for t in targets_1:
        if not re.match(multi_aa_rule.regex, t):
            targets.append(t)
        else:
            targets.extend(split_multiple_aa(t, multi_aa_rule.regex))
    result = old_coords_fix(coordinate_changes_dict[request.systematic_id], targets)
    return [OldCoordsFix.parse_obj(i) for i in result]


@ app.post("/histone_fix", response_model=list[AlleleFix])
async def fix_histone(request: HistoneFixRequest):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    targets_1 = request.targets.split(',')
    targets = list()
    for t in targets_1:
        if not re.match(multi_aa_rule.regex, t):
            targets.append(t)
        else:
            targets.extend(split_multiple_aa(t, multi_aa_rule.regex))
    print(request.dict() | {'targets': targets})
    result = apply_histone_fix(request.dict() | {'targets': ','.join(targets)}, genome, 'targets')
    return [AlleleFix.parse_obj({'values': result})] if result != '' else []


@app.get("/genome_region")
async def get_genome_region(systematic_id: str, format: str, upstream: int = 0, downstream: int = 0):

    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)

    if systematic_id not in genome:
        # For loci with multiple transcripts, we need to find the one that is the longest
        i = 1
        longest_transcript_systematic_id = None
        longest_transcript_length = 0
        while systematic_id + '.' + str(i) in genome:
            transcript_seq_record = extract_main_feature_and_strand(genome[systematic_id + '.' + str(i)], 0, 0)
            if len(transcript_seq_record) > longest_transcript_length:
                longest_transcript_length = len(transcript_seq_record)
                longest_transcript_systematic_id = systematic_id + '.' + str(i)
            i += 1
        if longest_transcript_systematic_id is None:
            raise HTTPException(404, 'Systematic id does not exist')
        else:
            systematic_id = longest_transcript_systematic_id

    seq_record, strand = extract_main_feature_and_strand(genome[systematic_id], downstream, upstream)
    seq_record.annotations['accession'] = systematic_id

    extension = 'gb'

    if format != 'genbank':
        extension = format

    with tempfile.NamedTemporaryFile('w', encoding='utf-8', delete=False) as fp:
        if strand == 1:
            SeqIO.write(seq_record, fp, format)
        else:
            rv = seq_record.reverse_complement()
            rv.annotations["molecule_type"] = "DNA"
            rv.annotations['accession'] = systematic_id
            SeqIO.write(rv, fp, format)

        return FileResponse(fp.name, background=BackgroundTask(lambda: os.remove(fp.name)), filename=f'{systematic_id}.{extension}')

@app.get('/residue_at_position')
async def get_residue_at_position(systematic_id: str, position: int, dna_or_protein: DNAorProtein):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    gene = genome[systematic_id]
    if dna_or_protein == 'dna':
        if 'CDS' in gene:
            seq = gene['CDS'].extract(gene['contig'])
        else:
            if len(gene) != 2:
                # Error, we cannot read this position
                raise ValueError('cannot read sequence, alternative splicing?')
            # The key is the one that is not 'contig'
            key = next(k for k in gene if k != 'contig')
            seq = gene[key].extract(gene['contig'])
        return seq.seq[position - 1]
    elif dna_or_protein == 'protein':
        if 'peptide' not in gene:
            raise ValueError('cannot read sequence, no peptide')
        return gene['peptide'][position - 1]