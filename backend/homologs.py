import re
import os
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from io import StringIO
from threading import Event

from Bio import Align
from Bio.Blast import NCBIXML

from backend.conservation import calculate_conservation, parse_aligned_fasta
from backend.functional_regions import fetch_functional_regions, fetch_uniprot_entry
from backend.mutation import clean_sequence


NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EBI_CLUSTALO_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_PEPTIDE_SEARCH_URL = "https://peptidesearch.uniprot.org/asyncrest/"
DEFAULT_CONTACT_EMAIL = os.getenv("MEWTATE_CONTACT_EMAIL", "brunotruzi@gmail.com")
NCBI_BLAST_TIMEOUT_SECONDS = 900
NCBI_BLAST_POLL_SECONDS = 30
UNIPROT_PEPTIDE_SEARCH_TIMEOUT_SECONDS = 20


def find_homologs_and_calculate_conservation(
    wildtype_sequence: str,
    max_homologs: int = 20,
    min_identity: float = 30.0,
    max_identity: float = 95.0,
    email: str | None = None,
    cancel_event: Event | None = None,
) -> dict:
    timings = {}
    total_start = time.perf_counter()
    check_cancelled(cancel_event)

    query_sequence = clean_sequence(wildtype_sequence)

    if not query_sequence:
        raise ValueError("Wildtype sequence cannot be empty.")

    contact_email = email or DEFAULT_CONTACT_EMAIL

    max_homologs = max(3, min(max_homologs, 50))

    fast_match_start = time.perf_counter()
    try:
        fast_match = find_uniprot_exact_match(query_sequence, cancel_event=cancel_event)
    except ValueError as error:
        fast_match = None
        timings["uniprot_exact_match_error"] = str(error)
        print(
            f"[deep-analysis] UniProt exact-match fast path failed: {error}. Falling back to NCBI BLAST.",
            flush=True,
        )

    timings["uniprot_exact_match_seconds"] = round(time.perf_counter() - fast_match_start, 2)

    if fast_match:
        uniprot_homolog_start = time.perf_counter()
        uniprot_homolog_result = build_uniprot_homolog_analysis(
            fast_match=fast_match,
            query_sequence=query_sequence,
            max_homologs=max_homologs,
            min_identity=min_identity,
            max_identity=max_identity,
            email=contact_email,
            cancel_event=cancel_event,
        )
        timings["uniprot_homolog_search_seconds"] = round(
            time.perf_counter() - uniprot_homolog_start,
            2,
        )

        if uniprot_homolog_result:
            timings.update(uniprot_homolog_result["timings"])
            timings["functional_regions"] = len(fast_match["functional_regions"])
            timings["total_seconds"] = round(time.perf_counter() - total_start, 2)
            print(
                "[deep-analysis timings] "
                + " ".join(f"{key}={value}" for key, value in timings.items()),
                flush=True,
            )

            return {
                **uniprot_homolog_result,
                "analysis_path": "uniprotkb_homolog_search",
                "functional_regions": fast_match["functional_regions"],
                "protein_match": fast_match["homologs"][0],
                "timings": timings,
                "filters": {
                    "database": "uniprotkb reviewed family/name search",
                    "identity_source": fast_match.get("database", "uniprotkb"),
                    "exact_sequence_match": True,
                    "max_homologs": max_homologs,
                    "min_identity": min_identity,
                    "max_identity": max_identity,
                    "remove_exact_duplicates": True,
                    "remove_fragments": True,
                    "prefer_reviewed": True,
                },
            }

        timings.update(fast_match["timings"])
        timings["total_seconds"] = round(time.perf_counter() - total_start, 2)
        print(
            "[deep-analysis timings] "
            + " ".join(f"{key}={value}" for key, value in timings.items()),
            flush=True,
        )

        return {
            **fast_match,
            "timings": timings,
            "filters": {
                "database": fast_match.get("database", "uniprotkb"),
                "exact_sequence_match": True,
                "prefer_reviewed": True,
                "homolog_search": "not enough UniProtKB homologs found",
            },
        }

    blast_start = time.perf_counter()
    blast_xml = run_ncbi_blastp(
        query_sequence,
        hitlist_size=50,
        cancel_event=cancel_event,
    )
    timings["ncbi_blast_seconds"] = round(time.perf_counter() - blast_start, 2)
    check_cancelled(cancel_event)

    parse_start = time.perf_counter()
    candidates = parse_blast_hits(
        blast_xml=blast_xml,
        query_length=len(query_sequence),
        min_identity=min_identity,
        max_identity=max_identity,
    )
    timings["parse_blast_hits_seconds"] = round(time.perf_counter() - parse_start, 2)
    timings["candidate_hits"] = len(candidates)

    fetch_start = time.perf_counter()
    homologs = fetch_and_filter_homolog_sequences(
        candidates=candidates,
        query_length=len(query_sequence),
        max_homologs=max_homologs,
        email=contact_email,
        cancel_event=cancel_event,
    )
    timings["fetch_filter_homologs_seconds"] = round(time.perf_counter() - fetch_start, 2)
    timings["filtered_homologs"] = len(homologs)
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

    msa_start = time.perf_counter()
    input_fasta = build_fasta(sequences_for_alignment)
    aligned_fasta = run_clustal_omega(
        input_fasta,
        email=contact_email,
        cancel_event=cancel_event,
    )
    timings["clustal_omega_seconds"] = round(time.perf_counter() - msa_start, 2)
    check_cancelled(cancel_event)

    conservation_start = time.perf_counter()
    aligned_sequences = parse_aligned_fasta(aligned_fasta)
    conservation = calculate_conservation(aligned_sequences)
    timings["conservation_seconds"] = round(time.perf_counter() - conservation_start, 2)

    uniprot_start = time.perf_counter()
    functional_regions = map_regions_to_query_coordinates(
        regions=fetch_functional_regions(homologs[0]["accession"]),
        top_homolog=homologs[0],
    )
    timings["uniprot_functional_regions_seconds"] = round(time.perf_counter() - uniprot_start, 2)
    timings["functional_regions"] = len(functional_regions)
    timings["total_seconds"] = round(time.perf_counter() - total_start, 2)

    print(
        "[deep-analysis timings] "
        + " ".join(f"{key}={value}" for key, value in timings.items()),
        flush=True,
    )

    return {
        "num_homologs": len(homologs),
        "homologs": homologs,
        "aligned_fasta": aligned_fasta,
        "conservation": conservation,
        "functional_regions": functional_regions,
        "timings": timings,
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


def find_uniprot_exact_match(
    query_sequence: str,
    cancel_event: Event | None = None,
) -> dict | None:
    exact_entry = find_uniprotkb_exact_sequence_match(
        query_sequence,
        cancel_event=cancel_event,
    )

    if exact_entry:
        return build_uniprot_exact_match_result(
            accession=exact_entry["accession"],
            entry=exact_entry["entry"],
            query_sequence=query_sequence,
            query_start=1,
            database="uniprotkb exact sequence search",
            cancel_event=cancel_event,
        )

    accessions = run_uniprot_peptide_search(query_sequence, cancel_event=cancel_event)

    for accession in accessions:
        check_cancelled(cancel_event)

        entry = fetch_uniprot_entry(accession)
        entry_sequence = entry.get("sequence", {}).get("value", "")
        query_offset = entry_sequence.find(query_sequence)

        if query_offset == -1:
            continue

        return build_uniprot_exact_match_result(
            accession=accession,
            entry=entry,
            query_sequence=query_sequence,
            query_start=query_offset + 1,
            database="uniprot peptide search",
            cancel_event=cancel_event,
        )

    return None


def find_uniprotkb_exact_sequence_match(
    query_sequence: str,
    cancel_event: Event | None = None,
) -> dict | None:
    check_cancelled(cancel_event)

    response = get_url(
        UNIPROT_SEARCH_URL,
        {
            "query": f"reviewed:true AND length:[{len(query_sequence)} TO {len(query_sequence)}]",
            "fields": "accession,protein_name,organism_name,length,sequence",
            "format": "json",
            "size": "500",
        },
    )
    data = json.loads(response)

    for entry in data.get("results", []):
        check_cancelled(cancel_event)

        entry_sequence = entry.get("sequence", {}).get("value", "")

        if entry_sequence == query_sequence:
            accession = entry.get("primaryAccession")

            if not accession:
                continue

            return {
                "accession": accession,
                "entry": fetch_uniprot_entry(accession),
            }

    return None


def build_uniprot_exact_match_result(
    accession: str,
    entry: dict,
    query_sequence: str,
    query_start: int,
    database: str,
    cancel_event: Event | None = None,
) -> dict:
    check_cancelled(cancel_event)

    query_end = query_start + len(query_sequence) - 1
    functional_regions = map_regions_to_query_subsequence(
        regions=fetch_functional_regions(accession),
        query_start=query_start,
        query_end=query_end,
    )
    homolog = build_uniprot_homolog_summary(
        accession=accession,
        entry=entry,
        query_sequence=query_sequence,
        query_start=query_start,
        query_end=query_end,
    )

    return {
        "analysis_path": "uniprot_exact_match",
        "database": database,
        "entry": entry,
        "num_homologs": 0,
        "homologs": [homolog],
        "aligned_fasta": "",
        "conservation": [],
        "functional_regions": functional_regions,
        "timings": {
            "functional_regions": len(functional_regions),
        },
    }


def build_uniprot_homolog_analysis(
    fast_match: dict,
    query_sequence: str,
    max_homologs: int,
    min_identity: float,
    max_identity: float,
    email: str | None,
    cancel_event: Event | None = None,
) -> dict | None:
    homologs = find_uniprotkb_homolog_candidates(
        entry=fast_match["entry"],
        query_sequence=query_sequence,
        max_homologs=max_homologs,
        min_identity=min_identity,
        max_identity=max_identity,
        cancel_event=cancel_event,
    )

    if len(homologs) < 2:
        return None

    sequences_for_alignment = [
        {
            "id": "wildtype_query",
            "description": "Input wildtype sequence",
            "sequence": query_sequence,
        },
        *homologs,
    ]

    msa_start = time.perf_counter()
    aligned_fasta = run_clustal_omega(
        build_fasta(sequences_for_alignment),
        email=email,
        cancel_event=cancel_event,
    )
    aligned_sequences = parse_aligned_fasta(aligned_fasta)
    conservation = calculate_conservation(aligned_sequences)

    return {
        "num_homologs": len(homologs),
        "homologs": homologs,
        "aligned_fasta": aligned_fasta,
        "conservation": conservation,
        "timings": {
            "filtered_uniprot_homologs": len(homologs),
            "clustal_omega_seconds": round(time.perf_counter() - msa_start, 2),
        },
    }


def find_uniprotkb_homolog_candidates(
    entry: dict,
    query_sequence: str,
    max_homologs: int,
    min_identity: float,
    max_identity: float,
    cancel_event: Event | None = None,
) -> list[dict]:
    search_queries = build_uniprot_homolog_queries(entry, query_sequence)

    if not search_queries:
        return []

    homologs = []
    seen_sequences = {query_sequence}
    seen_accessions = set()

    for search_query in search_queries:
        response = get_url(
            UNIPROT_SEARCH_URL,
            {
                "query": search_query,
                "fields": "accession,protein_name,organism_name,length,sequence",
                "format": "json",
                "size": "500",
            },
        )
        data = json.loads(response)

        for candidate in data.get("results", []):
            check_cancelled(cancel_event)

            accession = candidate.get("primaryAccession")
            sequence = candidate.get("sequence", {}).get("value", "")

            if not accession or not sequence:
                continue

            if accession in seen_accessions or sequence in seen_sequences:
                continue

            if "fragment" in extract_uniprot_protein_name(candidate).lower():
                continue

            identity = calculate_pairwise_identity(query_sequence, sequence)

            if identity < min_identity or identity > max_identity:
                continue

            seen_accessions.add(accession)
            seen_sequences.add(sequence)
            homologs.append({
                "accession": accession,
                "title": build_uniprot_candidate_title(candidate),
                "description": build_uniprot_candidate_title(candidate),
                "identity": identity,
                "coverage": round((min(len(query_sequence), len(sequence)) / len(query_sequence)) * 100, 2),
                "e_value": None,
                "bit_score": None,
                "id": safe_sequence_id(accession),
                "sequence": sequence,
                "length": len(sequence),
            })

            if len(homologs) >= max_homologs:
                break

        if len(homologs) >= max_homologs:
            break

    return homologs


def build_uniprot_homolog_queries(entry: dict, query_sequence: str) -> list[str]:
    length_margin = max(20, round(len(query_sequence) * 0.3))
    min_length = max(1, len(query_sequence) - length_margin)
    max_length = len(query_sequence) + length_margin
    length_filter = f"length:[{min_length} TO {max_length}]"
    queries = []

    for family_term in extract_uniprot_family_terms(entry):
        queries.append(f'reviewed:true AND "{family_term}" AND {length_filter}')

    protein_name = extract_uniprot_protein_name(entry)

    if protein_name:
        simplified_name = simplify_protein_name_for_search(protein_name)

        if simplified_name:
            queries.append(f'reviewed:true AND protein_name:"{simplified_name}" AND {length_filter}')

    return dedupe_preserving_order(queries)


def extract_uniprot_family_terms(entry: dict) -> list[str]:
    family_terms = []

    for comment in entry.get("comments", []):
        if comment.get("commentType") != "SIMILARITY":
            continue

        for text in comment.get("texts", []):
            value = text.get("value", "")
            match = re.search(r"belongs to the (.+?) family", value, flags=re.IGNORECASE)

            if not match:
                continue

            family_term = match.group(1).replace("-", " ").strip()

            if family_term:
                family_terms.append(family_term)

    return family_terms


def simplify_protein_name_for_search(protein_name: str) -> str:
    name = re.sub(r"\bprecursor\b", "", protein_name, flags=re.IGNORECASE)
    name = re.sub(r"\bisoform\b.*", "", name, flags=re.IGNORECASE)
    name = re.sub(r"\s+", " ", name).strip(" -")
    return name


def dedupe_preserving_order(values: list[str]) -> list[str]:
    seen = set()
    deduped = []

    for value in values:
        if value in seen:
            continue

        seen.add(value)
        deduped.append(value)

    return deduped


def calculate_pairwise_identity(query_sequence: str, subject_sequence: str) -> float:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(query_sequence, subject_sequence)[0]
    aligned_query = str(alignment[0])
    aligned_subject = str(alignment[1])
    comparable_positions = 0
    matches = 0

    for query_residue, subject_residue in zip(aligned_query, aligned_subject):
        if query_residue == "-" or subject_residue == "-":
            continue

        comparable_positions += 1

        if query_residue == subject_residue:
            matches += 1

    if comparable_positions == 0:
        return 0.0

    return round((matches / comparable_positions) * 100, 2)


def extract_uniprot_protein_name(entry: dict) -> str:
    return (
        entry.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value")
        or entry.get("proteinDescription", {})
        .get("submissionNames", [{}])[0]
        .get("fullName", {})
        .get("value")
        or ""
    )


def build_uniprot_candidate_title(entry: dict) -> str:
    protein_name = extract_uniprot_protein_name(entry)
    organism = entry.get("organism", {}).get("scientificName", "")
    accession = entry.get("primaryAccession", "")
    title = protein_name or accession

    if organism:
        title = f"{title} OS={organism}"

    return title


def run_uniprot_peptide_search(
    query_sequence: str,
    cancel_event: Event | None = None,
) -> list[str]:
    run_uniprot_peptide_search.last_location = None
    response = post_form(
        UNIPROT_PEPTIDE_SEARCH_URL,
        {
            "peps": query_sequence,
            "spOnly": "on",
        },
    )

    job_url = getattr(run_uniprot_peptide_search, "last_location", None)

    if not job_url:
        return []

    job_url = job_url.replace("http://", "https://")
    start_time = time.perf_counter()

    while time.perf_counter() - start_time < UNIPROT_PEPTIDE_SEARCH_TIMEOUT_SECONDS:
        check_cancelled(cancel_event)

        try:
            result = get_raw_url(job_url)
        except ValueError:
            sleep_with_cancel(5, cancel_event)
            continue

        return [
            accession.strip()
            for accession in result.split(",")
            if accession.strip()
        ]

    return []


def map_regions_to_query_subsequence(
    regions: list[dict],
    query_start: int,
    query_end: int,
) -> list[dict]:
    mapped_regions = []

    for region in regions:
        overlap_start = max(region["start"], query_start)
        overlap_end = min(region["end"], query_end)

        if overlap_start > overlap_end:
            continue

        mapped_regions.append({
            **region,
            "uniprot_start": region["start"],
            "uniprot_end": region["end"],
            "start": overlap_start - query_start + 1,
            "end": overlap_end - query_start + 1,
        })

    return mapped_regions


def build_uniprot_homolog_summary(
    accession: str,
    entry: dict,
    query_sequence: str,
    query_start: int,
    query_end: int,
) -> dict:
    protein_description = (
        entry.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value")
    )
    organism = entry.get("organism", {}).get("scientificName", "")
    title = protein_description or entry.get("uniProtkbId") or accession

    if organism:
        title = f"{title} OS={organism}"

    return {
        "accession": accession,
        "title": title,
        "description": title,
        "identity": 100.0,
        "coverage": 100.0,
        "e_value": 0,
        "bit_score": None,
        "id": safe_sequence_id(accession),
        "sequence": query_sequence,
        "length": len(query_sequence),
        "query_start": 1,
        "query_end": len(query_sequence),
        "subject_start": query_start,
        "subject_end": query_end,
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
    print(
        f"[deep-analysis] NCBI BLAST submitted RID={rid} estimated_wait={wait_seconds}s",
        flush=True,
    )
    sleep_with_cancel(min(wait_seconds, NCBI_BLAST_POLL_SECONDS), cancel_event)

    start_time = time.perf_counter()

    while time.perf_counter() - start_time < NCBI_BLAST_TIMEOUT_SECONDS:
        check_cancelled(cancel_event)

        status_response = get_url(
            NCBI_BLAST_URL,
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_OBJECT": "SearchInfo",
            },
        )
        elapsed_seconds = round(time.perf_counter() - start_time)
        status = extract_blast_value(status_response, "Status") or "UNKNOWN"

        print(
            f"[deep-analysis] NCBI BLAST status={status} elapsed={elapsed_seconds}s",
            flush=True,
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

        sleep_with_cancel(NCBI_BLAST_POLL_SECONDS, cancel_event)

    raise ValueError(
        f"NCBI BLAST search timed out after {NCBI_BLAST_TIMEOUT_SECONDS // 60} minutes."
    )


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
            "query_start": hsp.query_start,
            "query_end": hsp.query_end,
            "subject_start": hsp.sbjct_start,
            "subject_end": hsp.sbjct_end,
            "query_alignment": hsp.query,
            "subject_alignment": hsp.sbjct,
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


def map_regions_to_query_coordinates(
    regions: list[dict],
    top_homolog: dict,
) -> list[dict]:
    subject_to_query = build_subject_to_query_coordinate_map(top_homolog)
    mapped_regions = []

    for region in regions:
        mapped_start = subject_to_query.get(region["start"])
        mapped_end = subject_to_query.get(region["end"])

        if mapped_start is None and mapped_end is None:
            continue

        if mapped_start is None:
            mapped_start = mapped_end

        if mapped_end is None:
            mapped_end = mapped_start

        mapped_regions.append({
            **region,
            "uniprot_start": region["start"],
            "uniprot_end": region["end"],
            "start": min(mapped_start, mapped_end),
            "end": max(mapped_start, mapped_end),
        })

    return mapped_regions


def build_subject_to_query_coordinate_map(top_homolog: dict) -> dict[int, int]:
    query_position = top_homolog["query_start"] - 1
    subject_position = top_homolog["subject_start"] - 1
    subject_to_query = {}

    for query_residue, subject_residue in zip(
        top_homolog["query_alignment"],
        top_homolog["subject_alignment"],
    ):
        if query_residue != "-":
            query_position += 1

        if subject_residue != "-":
            subject_position += 1

        if query_residue != "-" and subject_residue != "-":
            subject_to_query[subject_position] = query_position

    return subject_to_query


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

    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            if url == UNIPROT_PEPTIDE_SEARCH_URL:
                run_uniprot_peptide_search.last_location = response.headers.get("Location")

            return response.read().decode("utf-8")
    except Exception as error:
        raise ValueError(f"Remote homolog search failed: {error}") from error


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


def get_raw_url(url: str) -> str:
    request = urllib.request.Request(
        url,
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
