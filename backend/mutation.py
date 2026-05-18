from backend.amino_acids import AMINO_ACIDS


def analyze_substitution(wildtype:str, mutant:str, position:int) -> dict:
    wildtype_props = AMINO_ACIDS[wildtype]
    mutant_props = AMINO_ACIDS[mutant]
    
    hydrophobicity_diff = (
        mutant_props["hydrophobicity"] - wildtype_props["hydrophobicity"]
    )
    
    severity_score = 0
    
    if wildtype_props["charge"] != mutant_props["charge"]:
        severity_score += 2
        
        
    if wildtype_props["polarity"] != mutant_props["polarity"]:
        severity_score += 1
        
    if abs(hydrophobicity_diff) >= 3:
        severity_score += 2
    elif abs(hydrophobicity_diff) >= 1.5:
        severity_score += 1
        
    if severity_score <= 1:
        severity = "low"
    elif severity_score <= 3:
        severity = "moderate"
    else:
        severity = "high"
        
        return{
            
            "type" : "substitution",
            "position": position,
            "wildtype" : wildtype,
            "wildtype_name" : wildtype_props["name"],
            "mutant" : mutant,
            "mutant_name" : mutant_props["name"],
            "mutation": f"{wildtype}{position}{mutant}",
            "charge_change": f"{wildtype_props['charge']} -> {mutant_props['charge']}",
            "hydrophobicity_change": round(hydrophobicity_diff,2),
            "severity": severity,
            "severity_score" : severity_score,
            
            }
    
    
    def compare_equal_length_proteins(
        wildtype_sequence:str,
        mutant_sequence:str,
        
    ) -> dict:
        wildtype_sequence = wildtype_sequence.upper().strip()
        mutant_sequence = mutant_sequence.upper().strip()
        
        
        if not wildtype_sequence:
            raise ValueError("Wildtype can't be empty")
        
        if not mutant_sequence:
            raise ValueError("Mutant can't be emtpy")
        
        if len(wildtype_sequence) != len(mutant_sequence):
            raise ValueError("Only equal lengths support at the moment")
        
        if set(wildtype_sequence) not in set(AMINO_ACIDS):
            raise ValueError("Wildtype contains invalid amino acids")
        
        if set(mutant_sequence) not in set(AMINO_ACIDS):
            raise ValueError("Mutant contains invalid amino acids")
        
        
