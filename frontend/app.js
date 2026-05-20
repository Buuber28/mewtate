const analyzeButton = document.getElementById("analyze-button");
const resultsSection = document.getElementById("results");

analyzeButton.addEventListener("click", async () => {
    const wildtypeSequence = document.getElementById("wildtype-sequence").value;
    const mutantSequence = document.getElementById("mutant-sequence").value;

    resultsSection.innerHTML = "<p>Analyzing...</p>";

    try {
        const response = await fetch("http://127.0.0.1:8000/compare-proteins", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({
                wildtype_sequence: wildtypeSequence,
                mutant_sequence: mutantSequence,
            }),
        });

        const data = await response.json();

        if (!response.ok) {
            resultsSection.innerHTML = `<p class="error">${data.detail}</p>`;
            return;
        }

        displayResults(data);
    } catch (error) {
        resultsSection.innerHTML = `<p class="error">Could not connect to backend.</p>`;
    }
});

function displayResults(data) {
    const substitutionsHtml = data.substitutions.map((substitution) => {
        return `
    <div class="substitution-card severity-${substitution.severity}">
        <div class="substitution-header">
            <strong>${substitution.mutation}</strong>
            <span>${substitution.severity}</span>
        </div>

        <p>${substitution.wildtype_name} → ${substitution.mutant_name}</p>

        <div class="substitution-details">
            <div><strong>Charge</strong><br>${substitution.charge_change}</div>
            <div><strong>Polarity</strong><br>${substitution.polarity_change}</div>
            <div><strong>Hydrophobicity</strong><br>${substitution.hydrophobicity_difference}</div>
            <div><strong>BLOSUM62</strong><br>${substitution.blosum62_score}</div>
        </div>
    </div>
`;
    }).join("");

    const sequenceVisualizationHtml = buildSequenceVisualization(data);

    resultsSection.innerHTML = `
        <h2>Results</h2>

        <p><strong>Sequence length:</strong> ${data.sequence_length}</p>
        <p><strong>Number of substitutions:</strong> ${data.num_substitutions}</p>
        <p><strong>Identity:</strong> ${data.identity_percent}%</p>
        <p><strong>Overall severity:</strong> ${data.overall_severity}</p>

        <h3>Sequence comparison</h3>
        ${sequenceVisualizationHtml}

        <h3>Detected substitutions</h3>
        <div class="substitution-list">
    ${substitutionsHtml || "<p>No substitutions detected.</p>"}
        </div>
    `;
}

function buildSequenceVisualization(data) {
    const wildtypeSequence = document.getElementById("wildtype-sequence").value.toUpperCase().trim();
    const mutantSequence = document.getElementById("mutant-sequence").value.toUpperCase().trim();

    const changedPositions = new Set(
        data.substitutions.map((substitution) => substitution.position)
    );

    let wildtypeTiles = "";
    let mutantTiles = "";

    for (let index = 0; index < wildtypeSequence.length; index++) {
        const position = index + 1;
        const isChanged = changedPositions.has(position);

        wildtypeTiles += `
            <span class="residue-tile ${isChanged ? "changed" : ""}" title="Position ${position}">
                ${wildtypeSequence[index]}
            </span>
        `;

        mutantTiles += `
            <span class="residue-tile ${isChanged ? "changed" : ""}" title="Position ${position}">
                ${mutantSequence[index]}
            </span>
        `;
    }

    return `
        <div class="sequence-visualization">
            <div class="sequence-row">
                <span class="sequence-label">WT</span>
                <div class="sequence-tiles">${wildtypeTiles}</div>
            </div>

            <div class="sequence-row">
                <span class="sequence-label">Mut</span>
                <div class="sequence-tiles">${mutantTiles}</div>
            </div>
        </div>
    `;
}