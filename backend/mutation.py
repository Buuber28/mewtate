from backend.amino_acids import AMINO_ACIDS
from backend.blosum import get_blosum62_score
from backend.alignment import (
    align_protein_sequences,
    calculate_identity,
)

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

    return edits

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
        return compare_equal_length_proteins(
            wildtype_sequence,
            mutant_sequence,
        )

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

    return {
        "sequence_length": len(wildtype_sequence),
        "identity_percent": identity,
        "alignment_score": alignment["score"],
        "aligned_wildtype": aligned_wildtype,
        "aligned_mutant": aligned_mutant,
        "edits": edits,
    }