from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware

from backend.conservation import parse_aligned_fasta, calculate_conservation
from backend.mutation import compare_proteins

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class SequenceComparisonRequest(BaseModel):
    wildtype_sequence: str
    mutant_sequence: str
    aligned_fasta: str | None = None

class ConservationRequest(BaseModel):
    aligned_fasta: str

@app.get("/")
def read_root():
    return {"message": "mewtate backend is running"}


@app.post("/compare-proteins")
def compare_proteins_endpoint(request: SequenceComparisonRequest):
    try:
        return compare_proteins(
    wildtype_sequence=request.wildtype_sequence,
    mutant_sequence=request.mutant_sequence,
    aligned_fasta=request.aligned_fasta,
    )
    except ValueError as error:
        raise HTTPException(status_code=400, detail=str(error))
    

@app.post("/conservation")
def conservation_endpoint(request: ConservationRequest):
    try:
        aligned_sequences = parse_aligned_fasta(request.aligned_fasta)
        conservation = calculate_conservation(aligned_sequences)

        return {
            "num_sequences": len(aligned_sequences),
            "alignment_length": len(aligned_sequences[0]),
            "conservation": conservation,
        }

    except ValueError as error:
        raise HTTPException(status_code=400, detail=str(error))