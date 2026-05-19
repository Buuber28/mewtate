from backend.amino_acids import AMINO_ACIDS


def analyze_substitution(wildtype: str, mutant: str, position: int) -> dict:
    wildtype_props = AMINO_ACIDS[wildtype]
    mutant_props = AMINO_ACIDS[mutant]

    hydrophobicity_difference = (
        mutant_props["hydrophobicity"] - wildtype_props["hydrophobicity"]
    )

    severity_score = 0

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
        "severity": severity,
        "severity_score": severity_score,
    }


def compare_equal_length_proteins(
    wildtype_sequence: str,
    mutant_sequence: str,
) -> dict:
    wildtype_sequence = wildtype_sequence.upper().strip()
    mutant_sequence = mutant_sequence.upper().strip()

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

    return {
        "sequence_length": len(wildtype_sequence),
        "num_substitutions": len(substitutions),
        "identity_percent": round(identity * 100, 2),
        "overall_severity": overall_severity,
        "substitutions": substitutions,
    }