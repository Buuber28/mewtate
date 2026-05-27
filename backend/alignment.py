from Bio import Align


def align_protein_sequences(wildtype_sequence: str, mutant_sequence: str) -> dict:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"

    aligner.match_score = 1  
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(wildtype_sequence, mutant_sequence)
    best_alignment = alignments[0]

    aligned_wildtype = str(best_alignment[0])
    aligned_mutant = str(best_alignment[1])

    return {
        "aligned_wildtype": aligned_wildtype,
        "aligned_mutant": aligned_mutant,
        "score": best_alignment.score,
    }


def calculate_identity(aligned_wildtype: str, aligned_mutant: str) -> float:
    comparable_positions = 0
    matches = 0

    for wt_residue, mut_residue in zip(aligned_wildtype, aligned_mutant):
        if wt_residue == "-" or mut_residue == "-":
            continue

        comparable_positions += 1

        if wt_residue == mut_residue:
            matches += 1

    if comparable_positions == 0:
        return 0.0

    return round((matches / comparable_positions) * 100, 2)