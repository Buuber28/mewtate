const analyzeButton = document.getElementById("analyze-button");
const resultsSection = document.getElementById("results");
const loadSequenceButton = document.getElementById("load-sequence-button");
const sequenceEditor = document.getElementById("sequence-editor");
const aminoAcidPalette = document.getElementById("amino-acid-palette");
const aminoAcidLegend = document.getElementById("amino-acid-legend");
const substituteModeButton = document.getElementById("substitute-mode-button");
const insertModeButton = document.getElementById("insert-mode-button");
const deleteModeButton = document.getElementById("delete-mode-button");
const findHomologsButton = document.getElementById("find-homologs-button");
const homologSearchStatus = document.getElementById("homolog-search-status");
const homologSearchProgress = document.getElementById("homolog-search-progress");
const proteinInfo = document.getElementById("protein-info");

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
let foundHomologAlignedFasta = "";
let foundHomologCount = 0;
let homologSearchController = null;



analyzeButton.addEventListener("click", async () => {
    const customAlignedFasta = document.getElementById("aligned-fasta").value.trim();
    const alignedFasta = customAlignedFasta || foundHomologAlignedFasta;
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

findHomologsButton.addEventListener("click", async () => {
    const wildtypeSequence = document.getElementById("wildtype-sequence").value;

    if (homologSearchController) {
        homologSearchController.abort();
    }

    homologSearchController = new AbortController();

    homologSearchStatus.textContent = "Searching NCBI BLAST and building alignment...";
    homologSearchStatus.classList.remove("error", "success");
    homologSearchProgress.classList.remove("hidden");
    findHomologsButton.disabled = true;

    try {
        const response = await fetch("http://127.0.0.1:8000/find-homologs", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({
                wildtype_sequence: wildtypeSequence,
                max_homologs: 20,
                min_identity: 30,
                max_identity: 95,
            }),
            signal: homologSearchController.signal,
        });

        const data = await response.json();

        if (!response.ok) {
            homologSearchStatus.textContent = data.detail;
            homologSearchStatus.classList.add("error");
            return;
        }

        foundHomologAlignedFasta = data.aligned_fasta;
        foundHomologCount = data.num_homologs;
        updateProteinInfo(data.homologs?.[0]);
        homologSearchStatus.textContent =
            `Homolog alignment ready with ${foundHomologCount} filtered homologs. Press Analyze mutation to use it.`;
        homologSearchStatus.classList.add("success");
    } catch (error) {
        if (error.name === "AbortError") {
            homologSearchStatus.textContent = "Homolog search cancelled.";
            return;
        }

        homologSearchStatus.textContent = "Could not connect to backend.";
        homologSearchStatus.classList.add("error");
    } finally {
        findHomologsButton.disabled = false;
        homologSearchProgress.classList.add("hidden");
        homologSearchController = null;
    }
});

window.addEventListener("pagehide", () => {
    if (homologSearchController) {
        homologSearchController.abort();
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
                        ${buildConservationDetails(edit)}
                    </div>
                </div>
            `;
        }

        if (edit.type === "insertion") {
            const insertionLabel = edit.length > 1
                ? `Insertion of ${edit.length} residues after position ${edit.position}`
                : `Insertion after position ${edit.position}`;
            const insertedText = edit.length > 1 ? "Inserted sequence" : "Inserted amino acid";

            return `
                <div class="substitution-card severity-${edit.severity}">
                    <div class="substitution-header">
                        <strong>${insertionLabel}</strong>
                        <span>${edit.severity}</span>
                    </div>
                    <p>${insertedText}: <strong>${edit.inserted}</strong></p>
                    ${buildIndelConservationDetails(edit)}
                </div>
            `;
        }

        if (edit.type === "deletion") {
            const deletionLabel = edit.length > 1
                ? `Deletion from positions ${edit.start_position}-${edit.end_position}`
                : `Deletion at position ${edit.position}`;
            const deletedText = edit.length > 1 ? "Deleted sequence" : "Deleted amino acid";

            return `
                <div class="substitution-card severity-${edit.severity}">
                    <div class="substitution-header">
                        <strong>${deletionLabel}</strong>
                        <span>${edit.severity}</span>
                    </div>
                    <p>${deletedText}: <strong>${edit.deleted}</strong></p>
                    ${buildIndelConservationDetails(edit)}
                </div>
            `;
        }

        return "";
    }).join("");

    const sequenceVisualizationHtml = buildSequenceVisualization(data);
    const conservationHeatmapHtml = buildConservationHeatmap(data);

    resultsSection.innerHTML = `
        <h2>Results</h2>

        <p><strong>Sequence length:</strong> ${data.sequence_length}</p>
        <p><strong>Identity:</strong> ${data.identity_percent}%</p>

        <h3>Sequence comparison</h3>
        ${sequenceVisualizationHtml}
        ${conservationHeatmapHtml}

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

function buildConservationDetails(edit) {
    if (!edit.conservation) {
        return `<div><strong>Conservation</strong><br>not provided</div>`;
    }

    return `
        <div><strong>Most common residue</strong><br>${edit.conservation.most_common_residue || "none"}</div>
        <div><strong>Conservation score</strong><br>${formatConservation(edit.conservation)}</div>
        <div><strong>Conservation label</strong><br>${edit.conservation.label}</div>
        <div><strong>Severity before conservation</strong><br>${edit.severity_before_conservation || edit.severity}</div>
        <div><strong>Severity after conservation</strong><br>${edit.severity_after_conservation || edit.severity}</div>
    `;
}

function buildIndelConservationDetails(edit) {
    if (!edit.conservation) {
        return "";
    }

    return `
        <div class="substitution-details">
            ${buildConservationDetails(edit)}
        </div>
    `;
}

function buildConservationHeatmap(data) {
    if (!data.conservation || data.conservation.length === 0) {
        return "";
    }

    const wildtypeSequence = cleanSequenceInput(
        document.getElementById("wildtype-sequence").value
    );

    const tiles = data.conservation
        .slice(0, wildtypeSequence.length)
        .map((item, index) => {
            const residue = wildtypeSequence[index] || item.most_common_residue || "-";
            const score = item.conservation_score;
            const heat = getConservationHeatColor(score);

            return `
                <span
                    class="conservation-tile"
                    style="background-color: ${heat.background}; border-color: ${heat.border}; color: ${heat.text};"
                    title="Position ${item.position}: ${score}% conserved, ${item.label}"
                >
                    ${residue}
                </span>
            `;
        })
        .join("");

    return `
        <h3>Conservation heatmap</h3>
        <div class="conservation-legend">
            <span>Variable</span>
            <span class="conservation-gradient"></span>
            <span>Highly conserved</span>
        </div>
        <div class="conservation-heatmap">${tiles}</div>
    `;
}

function getConservationHeatColor(score) {
    const normalizedScore = Math.max(0, Math.min(score, 100)) / 100;
    const lightness = 96 - (normalizedScore * 52);
    const saturation = 48 + (normalizedScore * 42);
    const text = normalizedScore >= 0.72 ? "white" : "#1f2937";

    return {
        background: `hsl(350, ${saturation}%, ${lightness}%)`,
        border: `hsl(350, ${Math.min(saturation + 6, 96)}%, ${Math.max(lightness - 12, 34)}%)`,
        text,
    };
}

function cleanSequenceInput(sequence) {
    return sequence
        .split("\n")
        .filter((line) => line.trim() && !line.trim().startsWith(">"))
        .join("")
        .replaceAll(" ", "")
        .toUpperCase();
}

function updateProteinInfo(topHomolog) {
    if (!topHomolog) {
        proteinInfo.innerHTML = `
            <p class="empty-state">
                No confident reviewed protein match was found.
            </p>
        `;
        return;
    }

    const confidence = getProteinMatchConfidence(topHomolog);

    proteinInfo.innerHTML = `
        <div class="protein-match">
            <span class="match-badge">${confidence}</span>
            <h3>${formatProteinTitle(topHomolog.description || topHomolog.title)}</h3>
            <dl>
                <div>
                    <dt>Accession</dt>
                    <dd>${topHomolog.accession}</dd>
                </div>
                <div>
                    <dt>Identity</dt>
                    <dd>${topHomolog.identity}%</dd>
                </div>
                <div>
                    <dt>Coverage</dt>
                    <dd>${topHomolog.coverage}%</dd>
                </div>
                <div>
                    <dt>Length</dt>
                    <dd>${topHomolog.length} aa</dd>
                </div>
            </dl>
        </div>
    `;
}

function getProteinMatchConfidence(topHomolog) {
    if (topHomolog.identity >= 98 && topHomolog.coverage >= 95) {
        return "Exact or near-exact match";
    }

    if (topHomolog.identity >= 80 && topHomolog.coverage >= 80) {
        return "Likely match";
    }

    return "Closest reviewed match";
}

function formatProteinTitle(title) {
    if (!title) {
        return "Reviewed Swiss-Prot protein";
    }

    const cleanedTitle = title
        .replace(/^>\s*/, "")
        .replace(/\s+OS=.*$/, "")
        .replace(/^sp\|[^|]+\|[^\s]+\s*/, "");

    const recommendedName = cleanedTitle.match(/RecName:\s*Full=([^;]+)/);

    if (recommendedName) {
        return recommendedName[1].trim();
    }

    return cleanedTitle
        .replace(/\s+Contains:.*$/, "")
        .replace(/\s+Flags:.*$/, "")
        .trim();
}

updateModeButtons();
