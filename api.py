from fastapi import FastAPI, HTTPException, Query
import json
from starlette.responses import RedirectResponse, PlainTextResponse, FileResponse
from pydantic import BaseModel
import pickle
from grammar import allowed_types_dict, composed_types_dict, aminoacid_grammar, nucleotide_grammar, disruption_grammar
from models import SyntaxRule, find_rule, AllowedTypes
from refinement_functions import check_allele_description, split_multiple_aa
from enum import Enum
from allele_fixes import multi_shift_fix, old_coords_fix, primer_mutagenesis as primer_mutagenesis_func
from common_autofix_functions import apply_histone_fix
from protein_modification_qc import check_func as check_modification_description
import re
from typing import Optional
from Bio import SeqIO
import tempfile
import os
from starlette.background import BackgroundTask
from genome_functions import extract_main_feature_and_strand, process_systematic_id, get_nt_at_gene_coord
from Bio.SeqRecord import SeqRecord
from transvar_functions import get_transvar_str_annotation, parse_transvar_string, TransvarAnnotation

syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]
syntax_rules_disruption = [SyntaxRule.parse_obj(r) for r in disruption_grammar]
multi_aa_rule = find_rule(syntax_rules_aminoacids, 'amino_acid_mutation', 'multiple_aa')
allowed_types = AllowedTypes(allowed_types=allowed_types_dict, composed_types=composed_types_dict)


class DNAorProtein(str, Enum):
    protein = 'protein'
    dna = 'dna'


class SequenceFileFormat(str, Enum):
    genbank = 'genbank'
    fasta = 'fasta'
    embl = 'embl'


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


# For the query field
systematic_id_description = 'Gene or transcript systematic id, if a gene that contains multiple transcripts is passed, the first transcript is used, `systematic_id.1`'
systematic_id_description_longest = 'Gene or transcript systematic id, if a gene that contains multiple transcripts is passed, the longest transcript is used'


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


def process_fix_targets(input_targets: list[str]):

    accepted_syntax = [r'[A-Z]?\d+[A-Z]?', r'\d+-\d+']

    targets = list()
    for t in input_targets:
        if re.fullmatch(multi_aa_rule.regex, t):
            # Format multi-substitution syntax
            targets.extend(split_multiple_aa(t, multi_aa_rule.regex))
        else:
            # Check all possible syntaxes
            if any(re.fullmatch(rule, t) for rule in accepted_syntax):
                targets.append(t)
            else:
                raise HTTPException(422, f'Invalid syntax for target {t}')

    return targets


def process_systematic_id_http_errors(systematic_id: str, genome: dict, when_several_transcripts: str) -> str:
    try:
        return process_systematic_id(systematic_id, genome, when_several_transcripts)
    except ValueError as e:
        raise HTTPException(404, str(e)) if 'Systematic id' in str(e) else HTTPException(400, str(e))


def extract_main_feature_and_strand_http_errors(gene: dict, downstream: int, upstream: int, get_utrs: bool = True) -> tuple[SeqRecord, int]:
    try:
        return extract_main_feature_and_strand(gene, downstream, upstream, get_utrs)
    except ValueError as e:
        raise HTTPException(400, str(e))


app = FastAPI()


@ app.get("/")
async def root():
    return RedirectResponse("/docs")


@ app.get("/check_allele", response_model=CheckAlleleDescriptionResponse)
async def check_allele(systematic_id: str = Query(example="SPBC359.03c", description=systematic_id_description), allele_description: str = Query(example="V123A,PLR-140-AAA,150-600"), allele_type: AlleleType = Query(example="partial_amino_acid_deletion")):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    systematic_id = process_systematic_id_http_errors(systematic_id, genome, 'first')
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
async def check_modification(systematic_id: str = Query(example="SPBC359.03c", description=systematic_id_description), sequence_position: str = Query(example="V123; V124,V125")):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    systematic_id = process_systematic_id_http_errors(systematic_id, genome, 'first')
    errors, change_sequence_position_to = check_modification_description({'systematic_id': systematic_id, 'sequence_position': sequence_position}, genome)
    needs_fixing = errors != '' or change_sequence_position_to != ''
    response_data = CheckModificationResponse(sequence_error=errors, change_sequence_position_to=change_sequence_position_to, needs_fixing=needs_fixing)
    response_data.user_friendly_fields = get_modification_user_friendly_fields(response_data)
    return response_data

@app.get("/primer")
async def primer_mutagenesis(systematic_id: str = Query(example="SPAPB1A10.09", description=systematic_id_description),
                             primer: str = Query(example="TTAGAGGTTATTAATTCCTAAGAAGAAGAAATTTTGG"),
                             max_mismatch: int = Query(description='The maximum amount of residues that are allowed to change'),
                             dna_or_protein: DNAorProtein = Query(), upstream: int = 0, downstream: int = 0
                             ):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    systematic_id = process_systematic_id_http_errors(systematic_id, genome, 'first')
    if dna_or_protein == 'protein':
        has_peptide = True
        gene = genome[systematic_id]
        seq = gene['CDS'].extract(gene['contig'])
    else:
        has_peptide = False
        seq, strand = extract_main_feature_and_strand(genome[systematic_id], downstream, upstream)
        if strand == -1:
            seq = seq.reverse_complement()
    primer = re.sub('\s+', '', primer.upper())
    print(primer)
    resp_str = primer_mutagenesis_func(seq.seq, primer, max_mismatch, has_peptide)
    if upstream != 0 and dna_or_protein == 'dna':
        substitutions = re.findall(r'[A-Z]-?\d+[A-Z]', resp_str)
        for substitution in substitutions:
            substitution_shifted = substitution.replace(substitution[1:-1],str(int(substitution[1:-1])-1))
            resp_str = resp_str.replace(substitution, substitution_shifted)

    return PlainTextResponse(resp_str)


# The same endpoint as above as a get endpoint
@ app.get("/multi_shift_fix", response_model=list[AlleleFix])
async def fix_with_multi_shift(systematic_id: str = Query(example="SPAPB1A10.09", description=systematic_id_description), targets: str = Query(example="S123,A124,N125"), dna_or_protein: DNAorProtein = DNAorProtein('protein')):

    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    systematic_id = process_systematic_id_http_errors(systematic_id, genome, 'first')

    gene = genome[systematic_id]
    targets = process_fix_targets(targets.split(','))

    if len([t for t in targets if t[0].isalpha()]) >= 3:
        if dna_or_protein == 'dna':
            seq, strand = extract_main_feature_and_strand(genome[systematic_id], 0, 0, get_utrs=False)
            if strand == -1:
                seq = seq.reverse_complement()

            result = multi_shift_fix(seq.seq, targets)
        else:
            result = multi_shift_fix(gene['peptide'], targets)
        return [AlleleFix.parse_obj({'values': i}) for i in result]
    else:
        return []


@ app.get("/old_coords_fix", response_model=list[OldCoordsFix])
async def fix_with_old_coords(systematic_id: str = Query(example="SPBC1706.01", description=systematic_id_description), targets: str = Query(example="P170A,V223A,F225A,AEY-171-LLL")):
    with open('data/coordinate_changes_dict.json') as ins:
        coordinate_changes_dict = json.load(ins)
    # We load the genome just to check if the systematic ID is valid
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    systematic_id = process_systematic_id_http_errors(systematic_id, genome, 'first')
    targets = process_fix_targets(targets.split(','))
    if systematic_id not in coordinate_changes_dict:
        return []
    result = old_coords_fix(coordinate_changes_dict[systematic_id], targets)
    return [OldCoordsFix.parse_obj(i) for i in result]


@ app.get("/histone_fix", response_model=list[AlleleFix])
async def fix_histone(systematic_id: str = Query(example="SPAC1834.04", description=systematic_id_description), targets: str = Query(example="ART-1-LLL,K9A,K14R,K14A")):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    systematic_id = process_systematic_id_http_errors(systematic_id, genome, 'first')
    targets = process_fix_targets(targets.split(','))
    result = apply_histone_fix({'systematic_id': systematic_id, 'targets': ','.join(targets)}, genome, 'targets')
    return [AlleleFix.parse_obj({'values': result})] if result != '' else []


@app.get("/genome_region")
async def get_genome_region(systematic_id: str = Query(example="SPAC1834.04", description=systematic_id_description_longest), format: SequenceFileFormat = Query(example="genbank"), upstream: int = 0, downstream: int = 0):

    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    systematic_id = process_systematic_id_http_errors(systematic_id, genome, 'longest')

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
async def get_residue_at_position(systematic_id: str = Query(example='SPAPB1A10.09', description=systematic_id_description), position: int = Query(example=1), dna_or_protein: DNAorProtein = Query(example='protein')):
    with open('data/genome.pickle', 'rb') as ins:
        genome = pickle.load(ins)
    if systematic_id not in genome:
        raise HTTPException(404, 'Systematic id does not exist')
    gene = genome[systematic_id]
    if dna_or_protein == 'dna':
        return PlainTextResponse(get_nt_at_gene_coord(position, gene, gene['contig']))
    elif dna_or_protein == 'protein':
        if 'peptide' not in gene:
            raise ValueError('cannot read sequence, no peptide')
        return PlainTextResponse(gene['peptide'][position - 1])


@ app.get("/ganno", summary='Variant described at the genome level (gDNA)', response_model=list[TransvarAnnotation])
async def ganno(variant_description: str = Query(example="II:g.178497T>A", description='Variant described at the genome level (gDNA)')) -> list[TransvarAnnotation]:
    try:
        return parse_transvar_string(get_transvar_str_annotation('ganno', variant_description))
    except Exception as e:
        raise HTTPException(400, str(e))


@ app.get("/canno", summary='Variant described at the coding DNA level (cDNA)', response_model=list[TransvarAnnotation])
async def canno(variant_description: str = Query(example="SPAC3F10.09:c.5A>T", description='Variant described at the coding DNA level (cDNA)')) -> list[TransvarAnnotation]:
    try:
        return parse_transvar_string(get_transvar_str_annotation('canno', variant_description))
    except Exception as e:
        raise HTTPException(400, str(e))


@ app.get("/panno", summary='Variant described at the protein level', response_model=list[TransvarAnnotation])
async def panno(variant_description: str = Query(example="SPBC1198.04c:p.N3A", description='Variant described at the protein level')) -> list[TransvarAnnotation]:
    try:
        return parse_transvar_string(get_transvar_str_annotation('panno', variant_description))
    except Exception as e:
        raise HTTPException(400, str(e))
