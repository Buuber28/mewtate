import asyncio
import contextlib
from threading import Event

from fastapi import FastAPI, HTTPException, Request
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware

from backend.conservation import parse_aligned_fasta, calculate_conservation
from backend.homologs import find_homologs_and_calculate_conservation
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
    functional_regions: list[dict] | None = None

class ConservationRequest(BaseModel):
    aligned_fasta: str

class HomologSearchRequest(BaseModel):
    wildtype_sequence: str
    max_homologs: int = 20
    min_identity: float = 30.0
    max_identity: float = 95.0
    email: str | None = None

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
            functional_regions=request.functional_regions,
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


@app.post("/find-homologs")
async def find_homologs_endpoint(
    http_request: Request,
    request: HomologSearchRequest,
):
    cancel_event = Event()

    async def watch_for_disconnect():
        while not cancel_event.is_set():
            if await http_request.is_disconnected():
                cancel_event.set()
                return

            await asyncio.sleep(1)

    disconnect_watcher = asyncio.create_task(watch_for_disconnect())

    try:
        return await asyncio.to_thread(
            find_homologs_and_calculate_conservation,
            wildtype_sequence=request.wildtype_sequence,
            max_homologs=request.max_homologs,
            min_identity=request.min_identity,
            max_identity=request.max_identity,
            email=request.email,
            cancel_event=cancel_event,
        )
    except ValueError as error:
        if cancel_event.is_set():
            raise HTTPException(status_code=499, detail="Homolog search was cancelled.")

        raise HTTPException(status_code=400, detail=str(error))
    finally:
        cancel_event.set()
        disconnect_watcher.cancel()

        with contextlib.suppress(asyncio.CancelledError):
            await disconnect_watcher
