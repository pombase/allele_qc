from fastapi import FastAPI
from starlette.responses import RedirectResponse, PlainTextResponse
from pydantic import BaseModel
import pickle
from grammar import allowed_types, aminoacid_grammar, nucleotide_grammar, disruption_grammar
from models import SyntaxRule
from refinement_functions import check_allele_description
from enum import Enum
from allele_fixes import multi_shift_fix, old_coords_fix
from typing import Optional

syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]
syntax_rules_disruption = [SyntaxRule.parse_obj(r) for r in disruption_grammar]


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


class CheckRequest(BaseModel):
    systematic_id: str
    allele_description: str
    allele_type: AlleleType

    class Config:
        schema_extra = {
            "example": {
                "systematic_id": "SPBC359.03c",
                "allele_description": "V123A,PLR-140-AAA,150-600",
                "allele_type": "amino_acid_deletion_and_mutation"}
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
                "targets": 'P170A,V223A,F225A'
            }
        }


class PrimerRequest(BaseModel):
    primer: str
    max_mismatch: int


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


def get_user_friendly_fields(obj: CheckAlleleDescriptionResponse) -> CheckAlleleDescriptionResponse:
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
            if 'position' in e:
                positions_dont_exist.append(e)
            else:
                residues_dont_match.append(e)
        out_str = ''
        if len(positions_dont_exist):
            out_str = out_str + f'The following sequence positions don\'t exist: {";".join(positions_dont_exist)}'
        if len(residues_dont_match):
            if len(positions_dont_exist):
                out_str = out_str + '\n'
            out_str = out_str + f'The following sequence positions don\'t contain the indicated residues: {";".join(residues_dont_match)}'
        new_obj.sequence_error = out_str
    if obj.change_type_to:
        new_obj.change_type_to = f'The indicated allele_type is not correct based on the existing mutations, it should be changed to {obj.change_type_to}'

    return new_obj


app = FastAPI()


@app.get("/")
async def root():
    return RedirectResponse("/docs")


@app.post("/check_allele", response_model=CheckAlleleDescriptionResponse)
async def check_allele(request: CheckRequest):
    with open('data/genome.pickle', 'rb') as ins:
        contig_genome = pickle.load(ins)

    response_data = CheckAlleleDescriptionResponse.parse_obj(
        check_allele_description(request.allele_description, syntax_rules_aminoacids, request.allele_type, allowed_types, contig_genome[request.systematic_id])
    )
    response_data.user_friendly_fields = get_user_friendly_fields(response_data)
    return response_data


@app.post("/primer")
async def primer_mutagenesis(request: CheckRequest):
    with open('data/genome.pickle', 'rb') as ins:
        contig_genome = pickle.load(ins)

    return check_allele_description(request.allele_description, syntax_rules_aminoacids, request.allele_type, allowed_types, contig_genome[request.systematic_id])


@app.post("/multi_shift")
async def multi_shift(request: MultiShiftRequest):
    with open('data/genome.pickle', 'rb') as ins:
        contig_genome = pickle.load(ins)
    gene = contig_genome[request.systematic_id]

    return PlainTextResponse(multi_shift_fix(gene['peptide'], request.targets.split(',')))


@app.post("/old_coords")
async def old_coords(request: OldCoordsFixRequest):

    return PlainTextResponse(old_coords_fix(request.systematic_id, request.targets.split(',')))
