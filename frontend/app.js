const quickAnalysisButton = document.getElementById("quick-analysis-button");
const deepAnalysisButton = document.getElementById("deep-analysis-button");
const deepAnalysisStatus = document.getElementById("deep-analysis-status");
const deepAnalysisProgress = document.getElementById("deep-analysis-progress");
const resultsSection = document.getElementById("results");
const recentSequences = document.getElementById("recent-sequences");
const sampleSequences = document.getElementById("sample-sequences");
const loadSequenceButton = document.getElementById("load-sequence-button");
const sequenceEditor = document.getElementById("sequence-editor");
const aminoAcidPalette = document.getElementById("amino-acid-palette");
const aminoAcidLegend = document.getElementById("amino-acid-legend");
const substituteModeButton = document.getElementById("substitute-mode-button");
const insertModeButton = document.getElementById("insert-mode-button");
const deleteModeButton = document.getElementById("delete-mode-button");
const proteinInfo = document.getElementById("protein-info");
const structurePanel = document.getElementById("structure-panel");
const structureViewerElement = document.getElementById("structure-viewer");
const structureStatus = document.getElementById("structure-status");
const structureExpandButton = document.getElementById("structure-expand-button");
const apiBaseUrl = "https://mewtate.onrender.com";

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

structureExpandButton?.addEventListener("click", () => {
    const isExpanded = structurePanel.classList.toggle("expanded");
    structureExpandButton.textContent = isExpanded ? "Collapse" : "Expand";
    refreshStructureViewer();
});

window.addEventListener("resize", () => {
    refreshStructureViewer();
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

const sampleProteinSequences = [
    {
        name: "Human insulin",
        detail: "INS precursor",
        sequence: "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
    },
    {
        name: "Hemoglobin beta",
        detail: "HBB mature chain",
        sequence: "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
    },
    {
        name: "Beta-lactamase SHV-1",
        detail: "P0AD63-like test",
        sequence: "MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW",
    },
];

let selectedResidueIndex = null;
let editableMutantSequence = "";

let selectedInsertPosition = null;

let mutationMode = "substitute";
let foundHomologAlignedFasta = "";
let foundHomologCount = 0;
let deepAnalysisController = null;
let foundFunctionalRegions = [];
let structureViewer = null;
let currentStructure = null;
const analysisCacheKey = "mewtate-analysis-cache";
const recentSequencesKey = "mewtate-recent-sequences";



quickAnalysisButton.addEventListener("click", async () => {
    await analyzeMutation({ useDeepAnalysis: false });
});

deepAnalysisButton.addEventListener("click", async () => {
    await analyzeMutation({ useDeepAnalysis: true });
});

async function analyzeMutation({ useDeepAnalysis }) {
    const customAlignedFasta = document.getElementById("aligned-fasta").value.trim();
    const wildtypeSequence = document.getElementById("wildtype-sequence").value;
    const mutantSequence = document.getElementById("mutant-sequence").value;
    let alignedFasta = customAlignedFasta;
    let functionalRegions = [];

    resultsSection.innerHTML = "<p>Analyzing...</p>";

    try {
        if (useDeepAnalysis) {
            const deepData = await runDeepAnalysisSetup(wildtypeSequence);
            alignedFasta = customAlignedFasta || deepData.aligned_fasta;
            functionalRegions = deepData.functional_regions || [];
        }

        const response = await fetch(`${apiBaseUrl}/compare-proteins`, {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({
                wildtype_sequence: wildtypeSequence,
                mutant_sequence: mutantSequence,
                aligned_fasta: alignedFasta || null,
                functional_regions: functionalRegions.length ? functionalRegions : null,
            }),
        });

        const data = await response.json();

        if (!response.ok) {
            resultsSection.innerHTML = `<p class="error">${data.detail}</p>`;
            return;
        }

        displayResults(data);
        updateStructureMutationHighlights(data.edits || data.substitutions || []);
    } catch (error) {
        if (error.name === "AbortError") {
            resultsSection.innerHTML = "<p>Deep analysis cancelled.</p>";
            return;
        }

        resultsSection.innerHTML = `<p class="error">Could not connect to backend.</p>`;
    } finally {
        deepAnalysisProgress.classList.add("hidden");
        quickAnalysisButton.disabled = false;
        deepAnalysisButton.disabled = false;
        deepAnalysisController = null;
    }
}

async function runDeepAnalysisSetup(wildtypeSequence) {
    const cleanedWildtype = cleanSequenceInput(wildtypeSequence);
    const cachedAnalysis = getCachedAnalysis(cleanedWildtype);

    if (cachedAnalysis) {
        applyDeepAnalysisData(cachedAnalysis);
        deepAnalysisStatus.textContent =
            `Reused cached deep analysis data for ${cachedAnalysis.display_name || "this sequence"} without rerunning BLAST/MSA. Running mutation analysis...`;
        deepAnalysisStatus.classList.add("success");
        return cachedAnalysis;
    }

    if (deepAnalysisController) {
        deepAnalysisController.abort();
    }

    deepAnalysisController = new AbortController();

    deepAnalysisStatus.textContent = "Finding homologs, building alignment, and fetching functional regions...";
    deepAnalysisStatus.classList.remove("error", "success");
    deepAnalysisProgress.classList.remove("hidden");
    quickAnalysisButton.disabled = true;
    deepAnalysisButton.disabled = true;

    const response = await fetch(`${apiBaseUrl}/find-homologs`, {
        method: "POST",
        headers: {
            "Content-Type": "application/json",
        },
        body: JSON.stringify({
            wildtype_sequence: cleanedWildtype,
            max_homologs: 10,
            min_identity: 30,
            max_identity: 95,
        }),
        signal: deepAnalysisController.signal,
    });

    const data = await response.json();

    if (!response.ok) {
        deepAnalysisStatus.textContent = data.detail;
        deepAnalysisStatus.classList.add("error");
        throw new Error(data.detail);
    }

    data.sequence = cleanedWildtype;
    data.display_name = getCachedDisplayName(data);
    cacheAnalysis(cleanedWildtype, data);
    applyDeepAnalysisData(data);
    renderRecentSequences();

    deepAnalysisStatus.textContent = getDeepAnalysisReadyMessage(data);
    deepAnalysisStatus.classList.add("success");

    return data;
}

function applyDeepAnalysisData(data) {
    foundHomologAlignedFasta = data.aligned_fasta || "";
    foundHomologCount = data.num_homologs || 0;
    foundFunctionalRegions = data.functional_regions || [];
    updateProteinInfo(data.protein_match || data.homologs?.[0], foundFunctionalRegions);
    loadProteinStructure(data.protein_match || data.homologs?.[0]);
}

window.addEventListener("pagehide", () => {
    if (deepAnalysisController) {
        deepAnalysisController.abort();
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
                        ${buildFunctionalRegionDetails(edit)}
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
                    ${buildIndelFunctionalRegionDetails(edit)}
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
                    ${buildIndelFunctionalRegionDetails(edit)}
                </div>
            `;
        }

        return "";
    }).join("");

    const sequenceVisualizationHtml = buildSequenceVisualization(data);
    const conservationHeatmapHtml = buildConservationHeatmap(data);
    const functionalRegionsHtml = buildFunctionalRegionsSummary(data);

    resultsSection.innerHTML = `
        <h2>Results</h2>

        <p><strong>Sequence length:</strong> ${data.sequence_length}</p>
        <p><strong>Identity:</strong> ${data.identity_percent}%</p>

        <h3>Sequence comparison</h3>
        ${sequenceVisualizationHtml}
        ${conservationHeatmapHtml}
        ${functionalRegionsHtml}

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

function buildFunctionalRegionDetails(edit) {
    if (!edit.functional_regions || edit.functional_regions.length === 0) {
        return `<div><strong>Functional region</strong><br>none detected</div>`;
    }

    const regionNames = edit.functional_regions
        .map((region) => `${region.description} (${region.start}-${region.end})`)
        .join("<br>");

    return `
        <div><strong>Functional region</strong><br>${regionNames}</div>
        <div><strong>Severity before functional region</strong><br>${edit.severity_before_functional_region || edit.severity}</div>
        <div><strong>Severity after functional region</strong><br>${edit.severity_after_functional_region || edit.severity}</div>
    `;
}

function buildIndelFunctionalRegionDetails(edit) {
    if (!edit.functional_regions || edit.functional_regions.length === 0) {
        return "";
    }

    return `
        <div class="substitution-details">
            ${buildFunctionalRegionDetails(edit)}
        </div>
    `;
}

function buildFunctionalRegionsSummary(data) {
    if (!data.functional_regions || data.functional_regions.length === 0) {
        return "";
    }

    const previewRegions = data.functional_regions.slice(0, 8).map((region) => `
        <span class="functional-region-chip">
            ${region.description} ${region.start}-${region.end}
        </span>
    `).join("");

    return `
        <h3>
            Functional regions
            <span class="info-tip" tabindex="0" data-tooltip="Functional regions come from UniProt feature annotations such as active sites, binding regions, domains, motifs, and transmembrane segments. Edits inside these regions can increase severity.">i</span>
        </h3>
        <div class="functional-region-list">
            ${previewRegions}
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
        <h3>
            Conservation heatmap
            <span class="info-tip" tabindex="0" data-tooltip="Darker red means the position is more conserved across homologs. A mutation in a highly conserved position is more likely to affect protein function.">i</span>
        </h3>
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

    return {
        background: `hsl(350, ${saturation}%, ${lightness}%)`,
        border: `hsl(350, ${Math.min(saturation + 6, 96)}%, ${Math.max(lightness - 12, 34)}%)`,
        text: "#111827",
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

function updateProteinInfo(topHomolog, functionalRegions = []) {
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
                <div>
                    <dt>Functional regions</dt>
                    <dd>${functionalRegions.length}</dd>
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

async function loadProteinStructure(proteinMatch) {
    if (!structureViewerElement || !structureStatus) {
        return;
    }

    if (!proteinMatch?.accession) {
        currentStructure = null;
        structureViewerElement.innerHTML = "";
        structureStatus.textContent = "Use Deep analysis to load an AlphaFold structure and map mutations onto it.";
        structureStatus.classList.remove("error", "success");
        return;
    }

    const accession = proteinMatch.accession;

    if (currentStructure?.accession === accession && structureViewer) {
        return;
    }

    const threeDmol = window.$3Dmol || window["3Dmol"];

    if (!threeDmol) {
        currentStructure = null;
        structureViewerElement.innerHTML = "";
        structureStatus.textContent = "3Dmol.js did not load, so the structure viewer is unavailable.";
        structureStatus.classList.add("error");
        structureStatus.classList.remove("success");
        return;
    }

    const structureUrls = [
        `https://alphafold.ebi.ac.uk/files/AF-${accession}-F1-model_v6.pdb`,
        `https://alphafold.ebi.ac.uk/files/AF-${accession}-F1-model_v5.pdb`,
        `https://alphafold.ebi.ac.uk/files/AF-${accession}-F1-model_v4.pdb`,
    ];
    structureStatus.textContent = `Loading AlphaFold model for ${accession}...`;
    structureStatus.classList.remove("error", "success");

    try {
        const structure = await fetchFirstAvailableStructure(structureUrls);
        structureViewerElement.innerHTML = "";
        structureViewer = threeDmol.createViewer(structureViewerElement, {
            backgroundColor: "#ffffff",
        });
        structureViewer.addModel(structure.pdbText, "pdb");
        structureViewer.setStyle({}, {
            cartoon: {
                color: "spectrum",
            },
        });
        structureViewer.zoomTo();
        structureViewer.render();
        refreshStructureViewer();

        currentStructure = {
            accession,
            structureUrl: structure.url,
            subjectStart: proteinMatch.subject_start || 1,
            queryLength: proteinMatch.length || null,
        };
        structureStatus.textContent =
            `AlphaFold model loaded for ${accession}. Run an analysis to highlight mutations.`;
        structureStatus.classList.add("success");
        structureStatus.classList.remove("error");
    } catch (error) {
        currentStructure = null;
        structureViewer = null;
        structureViewerElement.innerHTML = "";
        structureStatus.textContent =
            `No AlphaFold structure could be loaded for ${accession}.`;
        structureStatus.classList.add("error");
        structureStatus.classList.remove("success");
    }
}

async function fetchFirstAvailableStructure(urls) {
    for (const url of urls) {
        try {
            const response = await fetch(url);

            if (!response.ok) {
                continue;
            }

            return {
                url,
                pdbText: await response.text(),
            };
        } catch (error) {
            continue;
        }
    }

    throw new Error("AlphaFold model not found.");
}

function updateStructureMutationHighlights(edits) {
    if (!structureViewer || !currentStructure) {
        return;
    }

    const mutationPositions = getStructureMutationPositions(edits);

    structureViewer.setStyle({}, {
        cartoon: {
            color: "spectrum",
        },
    });

    if (mutationPositions.length) {
        structureViewer.addStyle(
            {
                resi: mutationPositions,
            },
            {
                stick: {
                    color: "#ef4444",
                    radius: 0.35,
                },
                sphere: {
                    color: "#ef4444",
                    radius: 0.75,
                },
            }
        );
        structureViewer.zoomTo({
            resi: mutationPositions,
        });
        structureStatus.textContent =
            `Highlighted ${mutationPositions.length} mutated structure position${mutationPositions.length === 1 ? "" : "s"} on ${currentStructure.accession}.`;
    } else {
        structureViewer.zoomTo();
        structureStatus.textContent =
            `AlphaFold model loaded for ${currentStructure.accession}. No residue-changing edits to highlight.`;
    }

    structureViewer.render();
    refreshStructureViewer();
}

function getStructureMutationPositions(edits) {
    const positions = new Set();

    edits.forEach((edit) => {
        if (edit.type === "insertion") {
            const insertionPosition = edit.position;

            if (!insertionPosition) {
                return;
            }

            const maxQueryPosition = currentStructure.queryLength;

            if (insertionPosition > 0) {
                positions.add(currentStructure.subjectStart + insertionPosition - 1);
            }

            if (!maxQueryPosition || insertionPosition < maxQueryPosition) {
                positions.add(currentStructure.subjectStart + insertionPosition);
            }

            return;
        }

        const start = edit.start_position || edit.position;
        const end = edit.end_position || edit.position;

        if (!start || !end) {
            return;
        }

        for (let position = start; position <= end; position++) {
            positions.add(currentStructure.subjectStart + position - 1);
        }
    });

    return [...positions];
}

function refreshStructureViewer() {
    if (!structureViewer) {
        return;
    }

    window.requestAnimationFrame(() => {
        if (typeof structureViewer.resize === "function") {
            structureViewer.resize();
        }

        structureViewer.render();
    });
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

function renderSequenceLibrary() {
    sampleSequences.innerHTML = sampleProteinSequences
        .map((sample, index) => buildSequenceButtonHtml(sample, index, "sample"))
        .join("");

    sampleSequences.querySelectorAll("button").forEach((button) => {
        button.addEventListener("click", () => {
            const sample = sampleProteinSequences[Number(button.dataset.index)];
            loadSequence(sample.sequence);
        });
    });

    renderRecentSequences();
}

function renderRecentSequences() {
    const recentItems = getRecentSequences();

    if (recentItems.length === 0) {
        recentSequences.innerHTML = `<p class="empty-state">Recent analyzed sequences will appear here.</p>`;
        return;
    }

    recentSequences.innerHTML = recentItems
        .map((item, index) => buildSequenceButtonHtml(item, index, "recent"))
        .join("");

    recentSequences.querySelectorAll("button").forEach((button) => {
        button.addEventListener("click", () => {
            const item = recentItems[Number(button.dataset.index)];
            loadSequence(item.sequence);

            const cachedAnalysis = getCachedAnalysis(item.sequence);

            if (cachedAnalysis) {
                applyDeepAnalysisData(cachedAnalysis);
                deepAnalysisStatus.textContent =
                    `Loaded cached deep analysis data for ${item.name}.`;
                deepAnalysisStatus.classList.add("success");
            }
        });
    });
}

function buildSequenceButtonHtml(item, index, source) {
    return `
        <button class="sequence-preset" data-source="${source}" data-index="${index}">
            <strong>${item.name}</strong>
            <span>${item.detail || `${item.sequence.length} aa`}</span>
        </button>
    `;
}

function loadSequence(sequence) {
    document.getElementById("wildtype-sequence").value = sequence;
    document.getElementById("mutant-sequence").value = sequence;
    aminoAcidLegend.classList.remove("hidden");
    editableMutantSequence = sequence;
    selectedResidueIndex = null;
    selectedInsertPosition = null;
    renderSequenceEditor();
    renderAminoAcidPalette();
}

function cacheAnalysis(sequence, data) {
    const cache = getAnalysisCache();
    cache[sequence] = data;
    sessionStorage.setItem(analysisCacheKey, JSON.stringify(cache));

    const recentItems = getRecentSequences()
        .filter((item) => item.sequence !== sequence);

    recentItems.unshift({
        name: data.display_name || "Analyzed protein",
        detail: getRecentSequenceDetail(sequence, data),
        sequence,
    });

    sessionStorage.setItem(
        recentSequencesKey,
        JSON.stringify(recentItems.slice(0, 6))
    );
}

function getCachedAnalysis(sequence) {
    return getAnalysisCache()[sequence] || null;
}

function getAnalysisCache() {
    try {
        return JSON.parse(sessionStorage.getItem(analysisCacheKey)) || {};
    } catch (error) {
        return {};
    }
}

function getRecentSequences() {
    try {
        return JSON.parse(sessionStorage.getItem(recentSequencesKey)) || [];
    } catch (error) {
        return [];
    }
}

function getCachedDisplayName(data) {
    const homolog = data.protein_match || data.homologs?.[0];

    if (!homolog) {
        return "Analyzed protein";
    }

    return formatProteinTitle(homolog.description || homolog.title);
}

function getRecentSequenceDetail(sequence, data) {
    if (data.analysis_path === "uniprotkb_homolog_search") {
        return `${sequence.length} aa · ${data.num_homologs || 0} UniProt homologs`;
    }

    if (data.analysis_path === "uniprot_exact_match") {
        return `${sequence.length} aa · UniProt exact match`;
    }

    return `${sequence.length} aa · ${data.num_homologs || 0} homologs`;
}

function formatTimingSummary(timings) {
    if (!timings) {
        return "";
    }

    if (timings.uniprot_homolog_search_seconds !== undefined) {
        return `Timing: UniProt match ${timings.uniprot_exact_match_seconds}s, homologs/MSA ${timings.uniprot_homolog_search_seconds}s.`;
    }

    if (timings.uniprot_exact_match_seconds !== undefined && timings.ncbi_blast_seconds === undefined) {
        return `Timing: UniProt exact match ${timings.uniprot_exact_match_seconds}s.`;
    }

    return `Timing: BLAST ${timings.ncbi_blast_seconds}s, fetch ${timings.fetch_filter_homologs_seconds}s, MSA ${timings.clustal_omega_seconds}s.`;
}

function getDeepAnalysisReadyMessage(data) {
    if (data.analysis_path === "uniprotkb_homolog_search") {
        return `Fast UniProt homolog analysis ready with ${data.num_homologs || 0} homologs. ${formatTimingSummary(data.timings)} Running mutation analysis...`;
    }

    if (data.analysis_path === "uniprot_exact_match") {
        return `Fast UniProt match found. ${formatTimingSummary(data.timings)} Running mutation analysis with protein identity and functional regions.`;
    }

    return `Deep analysis data ready with ${foundHomologCount} filtered homologs. ${formatTimingSummary(data.timings)} Running mutation analysis...`;
}

updateModeButtons();
renderSequenceLibrary();
