from backend.amino_acids import AMINO_ACIDS
from backend.blosum import get_blosum62_score
from backend.alignment import (
    align_protein_sequences,
    calculate_identity,
)
from backend.conservation import parse_aligned_fasta, calculate_conservation
from backend.functional_regions import apply_functional_regions_to_edits

def clean_sequence(sequence: str) -> str:
    lines = sequence.strip().splitlines()

    sequence_lines = [
        line.strip()
        for line in lines
        if line.strip() and not line.startswith(">")
    ]

    return "".join(sequence_lines).replace(" ", "").upper()

def analyze_substitution(wildtype: str, mutant: str, position: int) -> dict:
    wildtype_props = AMINO_ACIDS[wildtype]
    mutant_props = AMINO_ACIDS[mutant]

    hydrophobicity_difference = (
        mutant_props["hydrophobicity"] - wildtype_props["hydrophobicity"]
    )

    blosum62_score = get_blosum62_score(wildtype, mutant)

    severity_score = 0

    if blosum62_score <= -3:
        severity_score += 2
    elif blosum62_score < 0:
        severity_score += 1

    if wildtype_props["charge"] != mutant_props["charge"]:
        severity_score += 2

    if wildtype_props["polarity"] != mutant_props["polarity"]:
        severity_score += 1

    if abs(hydrophobicity_difference) >= 3:
        severity_score += 2
    elif abs(hydrophobicity_difference) >= 1.5:
        severity_score += 1

    if severity_score <= 1:
        severity = "low"
    elif severity_score <= 3:
        severity = "moderate"
    else:
        severity = "high"

    return {
        "position": position,
        "wildtype": wildtype,
        "wildtype_name": wildtype_props["name"],
        "mutant": mutant,
        "mutant_name": mutant_props["name"],
        "mutation": f"{wildtype}{position}{mutant}",
        "charge_change": f"{wildtype_props['charge']} → {mutant_props['charge']}",
        "polarity_change": f"{wildtype_props['polarity']} → {mutant_props['polarity']}",
        "hydrophobicity_difference": round(hydrophobicity_difference, 2),
        "blosum62_score": blosum62_score,
        "severity": severity,
        "severity_score": severity_score,
    }

def increase_severity(severity: str) -> str:
    severity_order = ["low", "moderate", "high"]

    if severity not in severity_order:
        return severity

    severity_index = severity_order.index(severity)
    next_index = min(severity_index + 1, len(severity_order) - 1)

    return severity_order[next_index]

def apply_conservation_to_edits(edits: list, conservation: list[dict]) -> None:
    conservation_by_position = {
        item["position"]: item
        for item in conservation
    }

    for edit in edits:
        affected_positions = get_edit_affected_positions(edit)
        conservation_items = [
            conservation_by_position[position]
            for position in affected_positions
            if position in conservation_by_position
        ]

        if not conservation_items:
            continue

        conservation_item = get_strongest_conservation_item(conservation_items)
        severity_before = edit["severity"]
        severity_after = severity_before

        if any(item["label"] == "highly conserved" for item in conservation_items):
            severity_after = increase_severity(severity_after)
        elif (
            any(item["label"] == "moderately conserved" for item in conservation_items)
            and severity_after == "low"
        ):
            severity_after = increase_severity(severity_after)

        edit["conservation"] = conservation_item
        edit["conservation_positions"] = conservation_items
        edit["severity_before_conservation"] = severity_before
        edit["severity_after_conservation"] = severity_after
        edit["severity"] = severity_after


def get_edit_affected_positions(edit: dict) -> list[int]:
    if edit.get("type") == "insertion":
        position = edit.get("position")

        if position is None:
            return []

        return [position, position + 1]

    start_position = edit.get("start_position", edit.get("position"))
    end_position = edit.get("end_position", edit.get("position"))

    if start_position is None or end_position is None:
        return []

    return list(range(start_position, end_position + 1))


def get_strongest_conservation_item(conservation_items: list[dict]) -> dict:
    label_rank = {
        "gap-only": 0,
        "variable": 1,
        "moderately conserved": 2,
        "highly conserved": 3,
    }

    return max(
        conservation_items,
        key=lambda item: (
            label_rank.get(item["label"], 0),
            item["conservation_score"],
        ),
    )

def detect_alignment_edits(
    aligned_wildtype: str,
    aligned_mutant: str,
) -> list:

    edits = []

    wildtype_position = 0

    for wt_residue, mut_residue in zip(aligned_wildtype, aligned_mutant):

        if wt_residue != "-":
            wildtype_position += 1

        # substitution
        if (
            wt_residue != "-"
            and mut_residue != "-"
            and wt_residue != mut_residue
        ):
            edits.append({
                "type": "substitution",
                **analyze_substitution(
                    wildtype=wt_residue,
                    mutant=mut_residue,
                    position=wildtype_position,
                )
            })

        # insertion
        elif wt_residue == "-" and mut_residue != "-":
            edits.append({
                "type": "insertion",
                "position": wildtype_position,
                "inserted": mut_residue,
                "severity": "moderate",
            })

        # deletion
        elif wt_residue != "-" and mut_residue == "-":
            edits.append({
                "type": "deletion",
                "position": wildtype_position,
                "deleted": wt_residue,
                "severity": "moderate",
            })

    return group_consecutive_indels(edits)

def group_consecutive_indels(edits: list) -> list:
    grouped_edits = []

    for edit in edits:
        if edit["type"] not in {"insertion", "deletion"}:
            grouped_edits.append(edit)
            continue

        previous_edit = grouped_edits[-1] if grouped_edits else None

        if edit["type"] == "insertion":
            can_group = (
                previous_edit
                and previous_edit["type"] == "insertion"
                and previous_edit["position"] == edit["position"]
            )

            if can_group:
                previous_edit["inserted"] += edit["inserted"]
                previous_edit["length"] += 1
            else:
                grouped_edits.append({
                    **edit,
                    "length": 1,
                })

            continue

        can_group = (
            previous_edit
            and previous_edit["type"] == "deletion"
            and previous_edit["end_position"] + 1 == edit["position"]
        )

        if can_group:
            previous_edit["deleted"] += edit["deleted"]
            previous_edit["end_position"] = edit["position"]
            previous_edit["length"] += 1
        else:
            grouped_edits.append({
                **edit,
                "start_position": edit["position"],
                "end_position": edit["position"],
                "length": 1,
            })

    return grouped_edits

def compare_equal_length_proteins(
    wildtype_sequence: str,
    mutant_sequence: str,
) -> dict:
    wildtype_sequence = clean_sequence(wildtype_sequence)
    mutant_sequence = clean_sequence(mutant_sequence)

    if not wildtype_sequence:
        raise ValueError("Wildtype sequence cannot be empty.")

    if not mutant_sequence:
        raise ValueError("Mutant sequence cannot be empty.")

    if len(wildtype_sequence) != len(mutant_sequence):
        raise ValueError(
            "Sequences must have the same length for this first version. "
        )

    invalid_wildtype = sorted(set(wildtype_sequence) - set(AMINO_ACIDS.keys()))
    invalid_mutant = sorted(set(mutant_sequence) - set(AMINO_ACIDS.keys()))

    if invalid_wildtype:
        raise ValueError(
            f"Wildtype sequence contains invalid amino acid symbols: {', '.join(invalid_wildtype)}"
        )

    if invalid_mutant:
        raise ValueError(
            f"Mutant sequence contains invalid amino acid symbols: {', '.join(invalid_mutant)}"
        )

    substitutions = []

    for index, (wildtype, mutant) in enumerate(
        zip(wildtype_sequence, mutant_sequence),
        start=1,
    ):
        if wildtype != mutant:
            substitutions.append(
                analyze_substitution(
                    wildtype=wildtype,
                    mutant=mutant,
                    position=index,
                )
            )

    if not substitutions:
        overall_severity = "none"
    else:
        max_score = max(change["severity_score"] for change in substitutions)

        if max_score <= 1:
            overall_severity = "low"
        elif max_score <= 3:
            overall_severity = "moderate"
        else:
            overall_severity = "high"

    identity = 1 - (len(substitutions) / len(wildtype_sequence))

    identity_percent = round(identity * 100, 2)

    if identity_percent < 70:
        raise ValueError(
            "Sequences are too different to interpret confidently."
        )

    return {
        "sequence_length": len(wildtype_sequence),
        "num_substitutions": len(substitutions),
        "identity_percent": round(identity * 100, 2),
        "overall_severity": overall_severity,
        "substitutions": substitutions,
    }

def compare_proteins(
    wildtype_sequence: str,
    mutant_sequence: str,
    aligned_fasta: str | None = None,
    functional_regions: list[dict] | None = None,
) -> dict:

    wildtype_sequence = clean_sequence(wildtype_sequence)
    mutant_sequence = clean_sequence(mutant_sequence)

    if not wildtype_sequence:
        raise ValueError("Wildtype sequence cannot be empty.")

    if not mutant_sequence:
        raise ValueError("Mutant sequence cannot be empty.")

    invalid_wildtype = sorted(set(wildtype_sequence) - set(AMINO_ACIDS.keys()))
    invalid_mutant = sorted(set(mutant_sequence) - set(AMINO_ACIDS.keys()))

    if invalid_wildtype:
        raise ValueError(
            f"Wildtype sequence contains invalid amino acid symbols: {', '.join(invalid_wildtype)}"
        )

    if invalid_mutant:
        raise ValueError(
            f"Mutant sequence contains invalid amino acid symbols: {', '.join(invalid_mutant)}"
        )

    # same-length fast path
    if len(wildtype_sequence) == len(mutant_sequence):
        result = compare_equal_length_proteins(
            wildtype_sequence,
            mutant_sequence,
        )

        edits = result["substitutions"]

        if aligned_fasta:
            aligned_sequences = parse_aligned_fasta(aligned_fasta)
            conservation = calculate_conservation(aligned_sequences)
            apply_conservation_to_edits(edits, conservation)
            result["conservation"] = conservation

        if functional_regions:
            apply_functional_regions_to_edits(edits, functional_regions)
            result["functional_regions"] = functional_regions

        result["edits"] = edits
        return result

    alignment = align_protein_sequences(
        wildtype_sequence,
        mutant_sequence,
    )

    aligned_wildtype = alignment["aligned_wildtype"]
    aligned_mutant = alignment["aligned_mutant"]

    identity = calculate_identity(
        aligned_wildtype,
        aligned_mutant,
    )

    if identity < 70:
        raise ValueError(
            "Sequences are too different to interpret confidently."
        )

    edits = detect_alignment_edits(
        aligned_wildtype,
        aligned_mutant,
    )

    if aligned_fasta:
        aligned_sequences = parse_aligned_fasta(aligned_fasta)
        conservation = calculate_conservation(aligned_sequences)
        apply_conservation_to_edits(edits, conservation)

    if functional_regions:
        apply_functional_regions_to_edits(edits, functional_regions)

    

    result = {
        "sequence_length": len(wildtype_sequence),
        "identity_percent": identity,
        "alignment_score": alignment["score"],
        "aligned_wildtype": aligned_wildtype,
        "aligned_mutant": aligned_mutant,
        "edits": edits,
    }

    if aligned_fasta:
        result["conservation"] = conservation

    if functional_regions:
        result["functional_regions"] = functional_regions

    return result
