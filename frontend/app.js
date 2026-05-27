const analyzeButton = document.getElementById("analyze-button");
const resultsSection = document.getElementById("results");
const loadSequenceButton = document.getElementById("load-sequence-button");
const sequenceEditor = document.getElementById("sequence-editor");
const aminoAcidPalette = document.getElementById("amino-acid-palette");
const aminoAcidLegend = document.getElementById("amino-acid-legend");
const substituteModeButton = document.getElementById("substitute-mode-button");
const insertModeButton = document.getElementById("insert-mode-button");
const deleteModeButton = document.getElementById("delete-mode-button");

substituteModeButton.addEventListener("click", () => {
    mutationMode = "substitute";
    resetSelections();
    updateModeButtons();
});

insertModeButton.addEventListener("click", () => {
    mutationMode = "insert";
    resetSelections();
    updateModeButtons();
});

deleteModeButton.addEventListener("click", () => {
    mutationMode = "delete";
    resetSelections();
    updateModeButtons();
});


const aminoAcids = [
    "A", "R", "N", "D", "C",
    "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P",
    "S", "T", "W", "Y", "V"
];

const aminoAcidClasses = {
    A: "nonpolar",
    V: "nonpolar",
    L: "nonpolar",
    I: "nonpolar",
    M: "nonpolar",
    F: "nonpolar",
    W: "nonpolar",
    P: "nonpolar",
    G: "nonpolar",

    S: "polar",
    T: "polar",
    N: "polar",
    Q: "polar",
    C: "polar",
    Y: "polar",

    K: "positive",
    R: "positive",
    H: "positive",

    D: "negative",
    E: "negative",
};

let selectedResidueIndex = null;
let editableMutantSequence = "";

let selectedInsertPosition = null;

let mutationMode = "substitute";



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

loadSequenceButton.addEventListener("click", () => {

   aminoAcidLegend.classList.remove("hidden");

    const wildtypeSequence = document
        .getElementById("wildtype-sequence")
        .value
        .toUpperCase()
        .trim();

    editableMutantSequence = wildtypeSequence;
    document.getElementById("mutant-sequence").value = editableMutantSequence;

    selectedResidueIndex = null;
    renderSequenceEditor();
    renderAminoAcidPalette();
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

function renderSequenceEditor() {
    sequenceEditor.innerHTML = "";

    if (mutationMode === "insert") {
        for (let index = 0; index <= editableMutantSequence.length; index++) {
            const insertSlot = document.createElement("button");
            insertSlot.textContent = "+";
            insertSlot.className = "insert-slot";

            if (index === selectedInsertPosition) {
                insertSlot.classList.add("selected-insert-slot");
            }

            insertSlot.title = `Insert at position ${index}`;

            insertSlot.addEventListener("click", () => {
                selectedInsertPosition = index;
                renderSequenceEditor();
            });

            sequenceEditor.appendChild(insertSlot);

            if (index < editableMutantSequence.length) {
                const residue = editableMutantSequence[index];

                const residueBlock = document.createElement("button");
                residueBlock.textContent = residue;
                residueBlock.className = `editable-residue aa-${aminoAcidClasses[residue]}`;
                residueBlock.disabled = true;

                sequenceEditor.appendChild(residueBlock);
            }
        }

        return;
    }

    for (let index = 0; index < editableMutantSequence.length; index++) {
        const residue = editableMutantSequence[index];

        const residueBlock = document.createElement("button");
        residueBlock.textContent = residue;
        residueBlock.className = `editable-residue aa-${aminoAcidClasses[residue]}`;
        residueBlock.title = `Position ${index + 1}`;

        if (mutationMode === "substitute" && index === selectedResidueIndex) {
            residueBlock.classList.add("selected");
        }

        residueBlock.addEventListener("click", () => {
            if (mutationMode === "substitute") {
                selectedResidueIndex = index;
            }

            if (mutationMode === "delete") {
                editableMutantSequence =
                    editableMutantSequence.slice(0, index) +
                    editableMutantSequence.slice(index + 1);

                document.getElementById("mutant-sequence").value = editableMutantSequence;
            }

            renderSequenceEditor();
        });

        sequenceEditor.appendChild(residueBlock);
    }
}

function renderAminoAcidPalette() {
    aminoAcidPalette.innerHTML = "";

    aminoAcids.forEach((aminoAcid) => {
        const aminoAcidBlock = document.createElement("button");
        aminoAcidBlock.textContent = aminoAcid;
        aminoAcidBlock.className = `amino-acid-block aa-${aminoAcidClasses[aminoAcid]}`;

        aminoAcidBlock.addEventListener("click", () => {
            if (mutationMode === "substitute") {
                if (selectedResidueIndex === null) {
                    return;
                }

                editableMutantSequence =
                    editableMutantSequence.slice(0, selectedResidueIndex) +
                    aminoAcid +
                    editableMutantSequence.slice(selectedResidueIndex + 1);
            }

            else if (mutationMode === "insert") {
                if (selectedInsertPosition === null) {
                    return;
                }

                editableMutantSequence =
                    editableMutantSequence.slice(0, selectedInsertPosition) +
                    aminoAcid +
                    editableMutantSequence.slice(selectedInsertPosition);
            }

            document.getElementById("mutant-sequence").value = editableMutantSequence;
            resetSelections();
        });

        aminoAcidPalette.appendChild(aminoAcidBlock);
    });
}

function resetSelections() {
    selectedResidueIndex = null;
    selectedInsertPosition = null;

    renderSequenceEditor();
}

function updateModeButtons() {
    substituteModeButton.classList.toggle("active-mode", mutationMode === "substitute");
    insertModeButton.classList.toggle("active-mode", mutationMode === "insert");
    deleteModeButton.classList.toggle("active-mode", mutationMode === "delete");
}