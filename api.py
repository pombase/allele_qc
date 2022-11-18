from fastapi import FastAPI
from starlette.responses import FileResponse
from pydantic import BaseModel
import pickle
from grammar import allowed_types, aminoacid_grammar, nucleotide_grammar, disruption_grammar
from models import SyntaxRule
from refinement_functions import check_allele_description


syntax_rules_aminoacids = [SyntaxRule.parse_obj(r) for r in aminoacid_grammar]
syntax_rules_nucleotides = [SyntaxRule.parse_obj(r) for r in nucleotide_grammar]
syntax_rules_disruption = [SyntaxRule.parse_obj(r) for r in disruption_grammar]


class MainRequest(BaseModel):
    systematic_id: str
    allele_description: str
    allele_type: str
    primer: str
    max_mismatch: int


app = FastAPI()


@app.get("/")
async def root():
    return FileResponse('static/index.html')


@app.post("/check_allele")
async def check_allele(request: MainRequest):
    with open('data/genome.pickle', 'rb') as ins:
        contig_genome = pickle.load(ins)

    return check_allele_description(request.allele_description, syntax_rules_aminoacids, request.allele_type, allowed_types, contig_genome[request.systematic_id])
