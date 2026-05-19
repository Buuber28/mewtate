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
            <li>
                <strong>${substitution.mutation}</strong>:
                ${substitution.wildtype_name} → ${substitution.mutant_name}
                <br>
                Charge: ${substitution.charge_change}
                <br>
                Polarity: ${substitution.polarity_change}
                <br>
                Hydrophobicity difference: ${substitution.hydrophobicity_difference}
                <br>
                Severity: ${substitution.severity}
            </li>
        `;
    }).join("");

    resultsSection.innerHTML = `
        <h2>Results</h2>

        <p><strong>Sequence length:</strong> ${data.sequence_length}</p>
        <p><strong>Number of substitutions:</strong> ${data.num_substitutions}</p>
        <p><strong>Identity:</strong> ${data.identity_percent}%</p>
        <p><strong>Overall severity:</strong> ${data.overall_severity}</p>

        <h3>Detected substitutions</h3>
        <ul>
            ${substitutionsHtml || "<li>No substitutions detected.</li>"}
        </ul>
    `;
}