from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from backend.mutation import compare_equal_length_proteins

app = FastAPI()


class SequenceComparisonRequest(BaseModel):
    wildtype_sequence: str
    mutant_sequence: str


@app.get("/")
def read_root():
    return {"message": "mewtate backend is running"}


@app.post("/compare-proteins")
def compare_proteins_endpoint(request: SequenceComparisonRequest):
    try:
        return compare_equal_length_proteins(
            wildtype_sequence=request.wildtype_sequence,
            mutant_sequence=request.mutant_sequence,
        )
    except ValueError as error:
        raise HTTPException(status_code=400, detail=str(error))