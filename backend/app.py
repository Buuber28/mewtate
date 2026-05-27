from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware

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


@app.get("/")
def read_root():
    return {"message": "mewtate backend is running"}


@app.post("/compare-proteins")
def compare_proteins_endpoint(request: SequenceComparisonRequest):
    try:
        return compare_proteins(
            wildtype_sequence=request.wildtype_sequence,
            mutant_sequence=request.mutant_sequence,
        )
    except ValueError as error:
        raise HTTPException(status_code=400, detail=str(error))