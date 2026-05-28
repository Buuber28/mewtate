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
    const alignedFasta = document.getElementById("aligned-fasta").value;
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
                aligned_fasta: alignedFasta || null,
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
    const edits = data.edits || data.substitutions || [];

    const editsHtml = edits.map((edit) => {
        if (edit.type === "substitution" || edit.mutation) {
            return `
                <div class="substitution-card severity-${edit.severity}">
                    <div class="substitution-header">
                        <strong>${edit.mutation}</strong>
                        <span>${edit.severity}</span>
                    </div>

                    <p>${edit.wildtype_name} → ${edit.mutant_name}</p>

                    <div class="substitution-details">
                        <div><strong>Charge</strong><br>${edit.charge_change}</div>
                        <div><strong>Polarity</strong><br>${edit.polarity_change}</div>
                        <div><strong>Hydrophobicity</strong><br>${edit.hydrophobicity_difference}</div>
                        <div><strong>BLOSUM62</strong><br>${edit.blosum62_score}</div>
                        <div><strong>Conservation</strong><br>${formatConservation(edit.conservation)}</div>
                    </div>
                </div>
            `;
        }

        if (edit.type === "insertion") {
            return `
                <div class="substitution-card severity-${edit.severity}">
                    <div class="substitution-header">
                        <strong>Insertion after position ${edit.position}</strong>
                        <span>${edit.severity}</span>
                    </div>
                    <p>Inserted amino acid: <strong>${edit.inserted}</strong></p>
                </div>
            `;
        }

        if (edit.type === "deletion") {
            return `
                <div class="substitution-card severity-${edit.severity}">
                    <div class="substitution-header">
                        <strong>Deletion at position ${edit.position}</strong>
                        <span>${edit.severity}</span>
                    </div>
                    <p>Deleted amino acid: <strong>${edit.deleted}</strong></p>
                </div>
            `;
        }

        return "";
    }).join("");

    const sequenceVisualizationHtml = buildSequenceVisualization(data);

    resultsSection.innerHTML = `
        <h2>Results</h2>

        <p><strong>Sequence length:</strong> ${data.sequence_length}</p>
        <p><strong>Identity:</strong> ${data.identity_percent}%</p>

        <h3>Sequence comparison</h3>
        ${sequenceVisualizationHtml}

        <h3>Detected edits</h3>
        <div class="substitution-list">
            ${editsHtml || "<p>No edits detected.</p>"}
        </div>
    `;
}

function buildSequenceVisualization(data) {
    const wildtypeSequence = data.aligned_wildtype
        || document.getElementById("wildtype-sequence").value.toUpperCase().trim();

    const mutantSequence = data.aligned_mutant
        || document.getElementById("mutant-sequence").value.toUpperCase().trim();

    const edits = data.edits || data.substitutions || [];

    const changedPositions = new Set(
        edits
            .filter((edit) => edit.type === "substitution" || edit.mutation)
            .map((edit) => edit.position)
    );

    const insertionPositions = new Set(
        edits
            .filter((edit) => edit.type === "insertion")
            .map((edit) => edit.position)
    );

    const deletionPositions = new Set(
        edits
            .filter((edit) => edit.type === "deletion")
            .map((edit) => edit.position)
    );

    let wildtypeTiles = "";
    let mutantTiles = "";
    let wildtypePosition = 0;

    for (let index = 0; index < wildtypeSequence.length; index++) {
        const wtResidue = wildtypeSequence[index];
        const mutResidue = mutantSequence[index];

        if (wtResidue !== "-") {
            wildtypePosition += 1;
        }

        const isSubstitution = changedPositions.has(wildtypePosition);
        const isInsertion = wtResidue === "-";
        const isDeletion = mutResidue === "-";

        wildtypeTiles += `
            <span class="residue-tile ${isSubstitution || isDeletion ? "changed" : ""} ${isInsertion ? "gap" : ""}" title="Alignment column ${index + 1}">
                ${wtResidue}
            </span>
        `;

        mutantTiles += `
            <span class="residue-tile ${isSubstitution || isInsertion ? "changed" : ""} ${isDeletion ? "gap" : ""}" title="Alignment column ${index + 1}">
                ${mutResidue}
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

function formatConservation(conservation) {
    if (!conservation) {
        return "not provided";
    }

    return `${conservation.conservation_score}%`;
}