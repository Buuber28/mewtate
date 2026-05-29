import re
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from io import StringIO
from threading import Event

from Bio.Blast import NCBIXML

from backend.conservation import calculate_conservation, parse_aligned_fasta
from backend.mutation import clean_sequence


NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EBI_CLUSTALO_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
DEFAULT_CONTACT_EMAIL = os.getenv("MEWTATE_CONTACT_EMAIL", "brunotruzi@gmail.com")


def find_homologs_and_calculate_conservation(
    wildtype_sequence: str,
    max_homologs: int = 20,
    min_identity: float = 30.0,
    max_identity: float = 95.0,
    email: str | None = None,
    cancel_event: Event | None = None,
) -> dict:
    check_cancelled(cancel_event)

    query_sequence = clean_sequence(wildtype_sequence)

    if not query_sequence:
        raise ValueError("Wildtype sequence cannot be empty.")

    contact_email = email or DEFAULT_CONTACT_EMAIL

    max_homologs = max(3, min(max_homologs, 50))

    blast_xml = run_ncbi_blastp(
        query_sequence,
        hitlist_size=50,
        cancel_event=cancel_event,
    )
    check_cancelled(cancel_event)

    candidates = parse_blast_hits(
        blast_xml=blast_xml,
        query_length=len(query_sequence),
        min_identity=min_identity,
        max_identity=max_identity,
    )

    homologs = fetch_and_filter_homolog_sequences(
        candidates=candidates,
        query_length=len(query_sequence),
        max_homologs=max_homologs,
        email=contact_email,
        cancel_event=cancel_event,
    )
    check_cancelled(cancel_event)

    if len(homologs) < 2:
        raise ValueError(
            "Not enough homologs passed the filters. Try a wider identity range."
        )

    sequences_for_alignment = [
        {
            "id": "wildtype_query",
            "description": "Input wildtype sequence",
            "sequence": query_sequence,
        },
        *homologs,
    ]

    input_fasta = build_fasta(sequences_for_alignment)
    aligned_fasta = run_clustal_omega(
        input_fasta,
        email=contact_email,
        cancel_event=cancel_event,
    )
    check_cancelled(cancel_event)

    aligned_sequences = parse_aligned_fasta(aligned_fasta)
    conservation = calculate_conservation(aligned_sequences)

    return {
        "num_homologs": len(homologs),
        "homologs": homologs,
        "aligned_fasta": aligned_fasta,
        "conservation": conservation,
        "filters": {
            "database": "swissprot",
            "max_homologs": max_homologs,
            "min_identity": min_identity,
            "max_identity": max_identity,
            "remove_exact_duplicates": True,
            "remove_fragments": True,
            "prefer_reviewed": True,
        },
    }


def run_ncbi_blastp(
    query_sequence: str,
    hitlist_size: int,
    cancel_event: Event | None = None,
) -> str:
    check_cancelled(cancel_event)

    submit_response = post_form(
        NCBI_BLAST_URL,
        {
            "CMD": "Put",
            "PROGRAM": "blastp",
            "DATABASE": "swissprot",
            "QUERY": query_sequence,
            "HITLIST_SIZE": str(hitlist_size),
            "EXPECT": "1e-5",
            "FILTER": "F",
        },
    )

    rid = extract_blast_value(submit_response, "RID")
    rtoe = extract_blast_value(submit_response, "RTOE")

    if not rid:
        raise ValueError("NCBI BLAST did not return a request ID.")

    wait_seconds = int(rtoe) if rtoe and rtoe.isdigit() else 10
    time.sleep(min(wait_seconds, 30))

    for _ in range(30):
        check_cancelled(cancel_event)

        status_response = get_url(
            NCBI_BLAST_URL,
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_OBJECT": "SearchInfo",
            },
        )

        if "Status=READY" in status_response:
            if "ThereAreHits=yes" not in status_response:
                raise ValueError("NCBI BLAST finished but found no homologs.")

            check_cancelled(cancel_event)

            return get_url(
                NCBI_BLAST_URL,
                {
                    "CMD": "Get",
                    "RID": rid,
                    "FORMAT_TYPE": "XML",
                },
            )

        if "Status=FAILED" in status_response:
            raise ValueError("NCBI BLAST search failed.")

        if "Status=UNKNOWN" in status_response:
            raise ValueError("NCBI BLAST request expired or was not found.")

        sleep_with_cancel(10, cancel_event)

    raise ValueError("NCBI BLAST search timed out.")


def parse_blast_hits(
    blast_xml: str,
    query_length: int,
    min_identity: float,
    max_identity: float,
) -> list[dict]:
    blast_record = NCBIXML.read(StringIO(blast_xml))
    candidates = []

    for alignment in blast_record.alignments:
        if not alignment.hsps:
            continue

        hsp = alignment.hsps[0]
        identity = round((hsp.identities / hsp.align_length) * 100, 2)
        coverage = round(((hsp.query_end - hsp.query_start + 1) / query_length) * 100, 2)

        if identity < min_identity or identity > max_identity:
            continue

        if coverage < 70:
            continue

        accession = normalize_accession(alignment.accession or alignment.hit_id)

        candidates.append({
            "accession": accession,
            "title": alignment.title,
            "identity": identity,
            "coverage": coverage,
            "e_value": hsp.expect,
            "bit_score": hsp.bits,
        })

    return candidates


def fetch_and_filter_homolog_sequences(
    candidates: list[dict],
    query_length: int,
    max_homologs: int,
    email: str | None,
    cancel_event: Event | None = None,
) -> list[dict]:
    homologs = []
    seen_sequences = set()

    for candidate in candidates:
        check_cancelled(cancel_event)

        if len(homologs) >= max_homologs:
            break

        fasta_text = fetch_protein_fasta(candidate["accession"], email=email)
        record = parse_single_fasta(fasta_text)

        if not record:
            continue

        sequence = record["sequence"]
        header = record["description"]

        if "fragment" in header.lower():
            continue

        length_ratio = len(sequence) / query_length

        if length_ratio < 0.7 or length_ratio > 1.3:
            continue

        if sequence in seen_sequences:
            continue

        seen_sequences.add(sequence)

        homologs.append({
            **candidate,
            "id": safe_sequence_id(candidate["accession"]),
            "description": header,
            "sequence": sequence,
            "length": len(sequence),
        })

    return homologs


def fetch_protein_fasta(accession: str, email: str | None) -> str:
    params = {
        "db": "protein",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "tool": "mewtate",
    }

    if email:
        params["email"] = email

    return get_url(NCBI_EFETCH_URL, params)


def run_clustal_omega(
    fasta_text: str,
    email: str | None,
    cancel_event: Event | None = None,
) -> str:
    check_cancelled(cancel_event)

    job_id = post_form(
        f"{EBI_CLUSTALO_URL}/run",
        {
            "email": email,
            "stype": "protein",
            "sequence": fasta_text,
            "outfmt": "fa",
        },
    ).strip()

    if not job_id:
        raise ValueError("Clustal Omega did not return a job ID.")

    for _ in range(60):
        check_cancelled(cancel_event)

        status = get_url(f"{EBI_CLUSTALO_URL}/status/{job_id}", {}).strip()

        if status == "FINISHED":
            result_type = find_clustal_fasta_result_type(job_id)
            result = get_url(f"{EBI_CLUSTALO_URL}/result/{job_id}/{result_type}", {})

            if result.lstrip().startswith(">"):
                return result

            return clustal_to_fasta(result)

        if status in {"ERROR", "FAILURE", "NOT_FOUND"}:
            raise ValueError(f"Clustal Omega job ended with status {status}.")

        sleep_with_cancel(5, cancel_event)

    raise ValueError("Clustal Omega alignment timed out.")


def check_cancelled(cancel_event: Event | None) -> None:
    if cancel_event and cancel_event.is_set():
        raise ValueError("Homolog search was cancelled.")


def sleep_with_cancel(seconds: int, cancel_event: Event | None) -> None:
    for _ in range(seconds):
        check_cancelled(cancel_event)
        time.sleep(1)


def find_clustal_fasta_result_type(job_id: str) -> str:
    result_types_xml = get_url(f"{EBI_CLUSTALO_URL}/resulttypes/{job_id}", {})
    root = ET.fromstring(result_types_xml)
    identifiers = [
        element.text
        for element in root.iter()
        if element.tag.endswith("identifier") and element.text
    ]

    for preferred_type in ("aln-fasta", "fasta", "fa", "out"):
        if preferred_type in identifiers:
            return preferred_type

    for identifier in identifiers:
        if "fasta" in identifier or identifier == "out":
            return identifier

    if "aln-clustal" in identifiers:
        return "aln-clustal"

    raise ValueError("Clustal Omega did not provide a usable alignment result.")


def clustal_to_fasta(clustal_text: str) -> str:
    sequences = {}

    for line in clustal_text.splitlines():
        if not line.strip() or line.startswith("CLUSTAL") or line.startswith("#"):
            continue

        parts = line.split()

        if len(parts) < 2:
            continue

        sequence_id, segment = parts[0], parts[1]

        if re.fullmatch(r"[*:.]+", sequence_id):
            continue

        sequences.setdefault(sequence_id, [])
        sequences[sequence_id].append(segment)

    if not sequences:
        raise ValueError("Could not parse Clustal Omega alignment output.")

    return build_fasta([
        {
            "id": sequence_id,
            "description": sequence_id,
            "sequence": "".join(segments),
        }
        for sequence_id, segments in sequences.items()
    ])


def build_fasta(records: list[dict]) -> str:
    fasta_parts = []

    for record in records:
        header = safe_sequence_id(record["id"])
        sequence = record["sequence"]
        lines = [sequence[index:index + 80] for index in range(0, len(sequence), 80)]
        fasta_parts.append(f">{header}\n" + "\n".join(lines))

    return "\n".join(fasta_parts)


def parse_single_fasta(fasta_text: str) -> dict | None:
    description = None
    sequence_lines = []

    for line in fasta_text.splitlines():
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            description = line[1:]
        else:
            sequence_lines.append(line)

    if not description or not sequence_lines:
        return None

    return {
        "description": description,
        "sequence": "".join(sequence_lines).upper(),
    }


def extract_blast_value(response_text: str, key: str) -> str | None:
    match = re.search(rf"{key}\s*=\s*(\S+)", response_text)
    return match.group(1) if match else None


def normalize_accession(accession: str) -> str:
    parts = accession.split("|")

    if len(parts) >= 2 and parts[1]:
        return parts[1]

    return accession


def safe_sequence_id(sequence_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", sequence_id)[:60]


def post_form(url: str, data: dict[str, str]) -> str:
    encoded_data = urllib.parse.urlencode(data).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=encoded_data,
        headers={"User-Agent": "mewtate/0.1"},
        method="POST",
    )

    return open_request(request)


def get_url(url: str, params: dict[str, str]) -> str:
    full_url = url

    if params:
        full_url = f"{url}?{urllib.parse.urlencode(params)}"

    request = urllib.request.Request(
        full_url,
        headers={"User-Agent": "mewtate/0.1"},
        method="GET",
    )

    return open_request(request)


def open_request(request: urllib.request.Request) -> str:
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return response.read().decode("utf-8")
    except Exception as error:
        raise ValueError(f"Remote homolog search failed: {error}") from error
