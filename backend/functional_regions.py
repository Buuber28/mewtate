import json
import urllib.parse
import urllib.request


UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb"

HIGH_IMPACT_FEATURES = {
    "active site",
    "binding site",
    "calcium binding",
    "dna binding",
    "metal binding",
    "nucleotide binding",
    "site",
    "zinc finger",
}

MODERATE_IMPACT_FEATURES = {
    "coiled coil",
    "compositionally biased region",
    "disulfide bond",
    "domain",
    "intramembrane",
    "modified residue",
    "motif",
    "region",
    "repeat",
    "signal",
    "topological domain",
    "transmembrane",
}

FUNCTIONAL_FEATURES = HIGH_IMPACT_FEATURES | MODERATE_IMPACT_FEATURES


def fetch_functional_regions(accession: str) -> list[dict]:
    entry = fetch_uniprot_entry(accession)
    regions = []

    for feature in entry.get("features", []):
        feature_type = feature.get("type")
        normalized_type = normalize_feature_type(feature_type)

        if normalized_type not in FUNCTIONAL_FEATURES:
            continue

        start = extract_feature_position(feature, "start")
        end = extract_feature_position(feature, "end")

        if start is None or end is None:
            continue

        regions.append({
            "type": feature_type,
            "normalized_type": normalized_type,
            "description": get_feature_description(feature, feature_type),
            "start": start,
            "end": end,
            "impact": get_feature_impact(normalized_type),
        })

    return regions


def apply_functional_regions_to_edits(edits: list, functional_regions: list[dict]) -> None:
    for edit in edits:
        affected_regions = [
            region
            for region in functional_regions
            if edit_overlaps_region(edit, region)
        ]

        if not affected_regions:
            continue

        severity_before = edit["severity"]
        severity_after = severity_before

        if any(region["impact"] == "high" for region in affected_regions):
            severity_after = increase_severity(severity_after)
        elif severity_after == "low":
            severity_after = increase_severity(severity_after)

        edit["functional_regions"] = affected_regions
        edit["severity_before_functional_region"] = severity_before
        edit["severity_after_functional_region"] = severity_after
        edit["severity"] = severity_after


def edit_overlaps_region(edit: dict, region: dict) -> bool:
    if edit.get("type") == "insertion":
        insertion_position = edit.get("position")

        if insertion_position is None:
            return False

        affected_positions = {insertion_position, insertion_position + 1}

        return any(
            region["start"] <= position <= region["end"]
            for position in affected_positions
        )

    edit_start = edit.get("start_position", edit.get("position"))
    edit_end = edit.get("end_position", edit.get("position"))

    if edit_start is None or edit_end is None:
        return False

    if region.get("normalized_type") == "disulfide bond":
        return (
            edit_start <= region["start"] <= edit_end
            or edit_start <= region["end"] <= edit_end
        )

    return edit_start <= region["end"] and edit_end >= region["start"]


def fetch_uniprot_entry(accession: str) -> dict:
    url = f"{UNIPROT_ENTRY_URL}/{urllib.parse.quote(accession)}.json"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "mewtate/0.1"},
        method="GET",
    )

    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return json.loads(response.read().decode("utf-8"))
    except Exception as error:
        raise ValueError(f"Could not fetch UniProt functional regions: {error}") from error


def extract_feature_position(feature: dict, key: str) -> int | None:
    position = feature.get("location", {}).get(key, {})
    value = position.get("value")

    if value is None:
        return None

    return int(value)


def get_feature_impact(feature_type: str) -> str:
    if feature_type in HIGH_IMPACT_FEATURES:
        return "high"

    return "moderate"


def normalize_feature_type(feature_type: str | None) -> str:
    if not feature_type:
        return ""

    return feature_type.replace("_", " ").strip().lower()


def get_feature_description(feature: dict, feature_type: str | None) -> str:
    description = feature.get("description")

    if description:
        return description

    ligand = feature.get("ligand", {}).get("name")

    if ligand:
        return f"{feature_type}: {ligand}"

    return feature_type or "Functional feature"


def increase_severity(severity: str) -> str:
    severity_order = ["low", "moderate", "high"]

    if severity not in severity_order:
        return severity

    severity_index = severity_order.index(severity)
    next_index = min(severity_index + 1, len(severity_order) - 1)

    return severity_order[next_index]
