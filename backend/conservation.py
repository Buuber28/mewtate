from collections import Counter


def parse_aligned_fasta(fasta_text: str) -> list[str]:
    sequences = []
    current_sequence = []

    for line in fasta_text.strip().splitlines():
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if current_sequence:
                sequences.append("".join(current_sequence).upper())
                current_sequence = []
        else:
            current_sequence.append(line.upper())

    if current_sequence:
        sequences.append("".join(current_sequence).upper())

    if not sequences:
        raise ValueError("No aligned sequences found.")

    alignment_length = len(sequences[0])

    if any(len(sequence) != alignment_length for sequence in sequences):
        raise ValueError("All aligned sequences must have the same length.")

    return sequences


def calculate_conservation(aligned_sequences: list[str]) -> list[dict]:
    conservation = []

    alignment_length = len(aligned_sequences[0])

    for index in range(alignment_length):
        column = [
            sequence[index]
            for sequence in aligned_sequences
            if sequence[index] != "-"
        ]

        if not column:
            conservation.append({
                "position": index + 1,
                "most_common_residue": None,
                "conservation_score": 0.0,
                "label": "gap-only",
            })
            continue

        counts = Counter(column)
        most_common_residue, count = counts.most_common(1)[0]
        score = count / len(column)

        if score >= 0.9:
            label = "highly conserved"
        elif score >= 0.7:
            label = "moderately conserved"
        else:
            label = "variable"

        conservation.append({
            "position": index + 1,
            "most_common_residue": most_common_residue,
            "conservation_score": round(score * 100, 2),
            "label": label,
        })

    return conservation