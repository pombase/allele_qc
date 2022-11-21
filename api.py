from fastapi import FastAPI
from starlette.responses import RedirectResponse, PlainTextResponse
from pydantic import BaseModel
import pickle
from grammar import allowed_types, aminoacid_grammar, nucleotide_grammar, disruption_grammar
from models import SyntaxRule
from refinement_functions import check_allele_description
from enum import Enum
from allele_fixes import multi_shift_fix, old_coords_fix

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


app = FastAPI()


@app.get("/")
async def root():
    return RedirectResponse("/docs")


@app.post("/check_allele")
async def check_allele(request: CheckRequest):
    with open('data/genome.pickle', 'rb') as ins:
        contig_genome = pickle.load(ins)

    return check_allele_description(request.allele_description, syntax_rules_aminoacids, request.allele_type, allowed_types, contig_genome[request.systematic_id])


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
