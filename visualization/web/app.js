import {chapters, initialState, percentage, reconcile, rows, suppressed} from "./state.js";
import {MolecularStage} from "./viewer.js";

const $ = selector => document.querySelector(selector);
const escape = value => String(value ?? "").replace(/[&<>"']/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));
const title = value => value.replaceAll("_", " ");
const shortSource = {boltzgen_consensus: "BoltzGen consensus", p2rank_geometry: "Pocket geometry", known_ligand: "Known-ligand reference", hybrid: "Hybrid signature"};
let bundle, manifest, state, stage;

function sourceBadge() {
  return `<span class="tag ${manifest.data_origin === "computed" ? "" : "warning"}">${escape(title(manifest.data_origin).toUpperCase())}</span>`;
}

function candidate() {
  return manifest.candidates.find(c => c.id === state.candidate);
}

function shell() {
  $("#app").innerHTML = `
    <aside class="sidebar"><div><div class="brand"><span class="brand-mark">A</span>AutoRepurpose</div><div class="edition">Research atlas / v0.2</div></div>
      <div class="nav-label">The research journey</div><nav aria-label="Presentation chapters">${chapters.map((chapter, i) => `<button class="nav-button" data-chapter="${i}"><span class="nav-number">0${i + 1}</span><span><span class="nav-title">${chapter[0]}</span><span class="nav-subtitle">${chapter[1]}</span></span></button>`).join("")}</nav>
      <div class="sidebar-bottom"><strong>Local research snapshot</strong>No live jobs. No API credentials.<br>Evidence stays with the result.</div></aside>
    <div class="main-wrap"><header class="topbar"><div class="breadcrumb">Research atlas &nbsp; / &nbsp; <b>${escape(manifest.target.name)}</b></div><div class="toolbar">${sourceBadge()}<a class="button" href="report.html" target="_blank" rel="noopener">Open report</a><a class="button" href="bundle.json" download>Export data</a></div></header>
    <main id="main" class="content" tabindex="-1"><div id="origin-banner"></div><div id="hero"></div>
      <div class="stage-layout"><section class="panel" aria-label="Molecular comparison"><div class="panel-heading"><h2>Interface explorer</h2><div class="segmented" role="group" aria-label="Molecular view">${[["context", "Context"], ["site", "Site"], ["frequency", "Frequency"], ["coverage", "Coverage"]].map(([id, text]) => `<button data-view="${id}">${text}</button>`).join("")}</div></div>
        <div class="viewer-grid"><div class="molecular-column"><div class="molecular-label">Target-side reference<b id="reference-label"></b></div><div id="viewer-reference" class="viewer" role="img" aria-label="Reference target structure; equivalent residue data below"></div></div><div class="molecular-column"><div class="molecular-label">Selected molecule<b id="candidate-label"></b></div><div id="viewer-candidate" class="viewer" role="img" aria-label="Selected candidate structure; equivalent residue data below"></div></div></div>
        <div id="legend" class="legend"></div><div id="viewer-status" class="viewer-status" role="status"></div></section>
      <aside id="inspector" class="panel inspector" aria-label="Selected result explanation"></aside></div>
      <div id="chapter-content" class="section"></div><div id="residue-selection" class="subtle" role="status"></div>
      <div class="bottom-nav"><button id="previous" class="button">Previous chapter</button><button id="next" class="button primary">Next chapter</button></div>
      <footer class="footer"><span>${escape(bundle.disclaimer)}</span><span>Run ${escape(manifest.run_id)} / schema ${escape(manifest.schema_version)}</span></footer>
    </main></div>`;
  if (manifest.data_origin !== "computed") {
    $("#origin-banner").innerHTML = `<div class="origin-banner"><strong>${escape(title(manifest.data_origin).toUpperCase())}</strong><span>${manifest.data_origin === "synthetic_fixture" ? "Synthetic coordinates and scores for interface testing only. No scientific result is represented." : "Reference structure example only. This is not an evaluated repurposing run."}</span></div>`;
  }
  stage = new MolecularStage($("#viewer-reference"), $("#viewer-candidate"), $("#viewer-status"), selectResidue);
  document.addEventListener("click", event => {
    const element = event.target.closest("button");
    if (!element) return;
    if (element.dataset.chapter !== undefined) changeChapter(Number(element.dataset.chapter));
    if (element.dataset.view) {state.view = element.dataset.view; render();}
    if (element.dataset.candidate) {state.candidate = element.dataset.candidate; state.view = "coverage"; render();}
    if (element.dataset.residue) selectResidue(element.dataset.residue);
  });
  $("#previous").addEventListener("click", () => changeChapter(Math.max(0, state.chapter - 1)));
  $("#next").addEventListener("click", () => changeChapter(Math.min(5, state.chapter + 1)));
}

function changeChapter(index) {
  state.chapter = index;
  state.view = ["context", "site", "frequency", "coverage", "coverage", "frequency"][index];
  render();
}

function selectResidue(key) {
  state.residue = state.residue === key ? null : key;
  render();
}

function inspector() {
  const c = candidate();
  const sig = manifest.signature;
  const core = sig.residues.filter(r => r.core);
  if (state.chapter < 3) {
    return `<div class="eyebrow">${state.chapter === 0 ? "The starting point" : "The interface"}</div><h2>${escape(manifest.target.name)}</h2><p>${escape(manifest.target.site_kind)} · ${escape(manifest.target.id)}</p><div class="score">${sig.status === "available" ? core.length : "—"}</div><div class="score-label">${sig.status === "available" ? "consensus-defined core residues" : "Signature not available"}</div><hr><div class="detail-row"><span>Signature source</span><b>${escape(shortSource[sig.source])}</b></div><div class="detail-row"><span>Core threshold</span><b>${Math.round(sig.threshold * 100)}%</b></div><div class="detail-row"><span>Addressability</span><b>${escape(title(sig.addressability.status))}</b></div><hr><p>${escape(sig.addressability.reason)}</p><small>Consensus contact frequency is not evidence of therapeutic efficacy.</small>`;
  }
  if (!c) return `<div class="eyebrow">Selected hypothesis</div><h2>No candidate selected</h2><p>No records match these filters. Missing or rejected results are never replaced with invented hits.</p>`;
  const coverage = percentage(c.coverage);
  return `<div class="eyebrow">Selected hypothesis</div><h2>${escape(c.name)}</h2><span class="tag ${c.confidence.gate === "passed" ? "" : "warning"}">${escape(title(c.novelty))}</span><div class="score">${escape(coverage)}</div><div class="score-label">${c.coverage.status === "available" ? `${c.engaged.length} of ${core.length} core residues engaged` : escape(c.coverage.reason)}</div><hr><div class="detail-row"><span>Confidence gate</span><b>${escape(title(c.confidence.gate))}</b></div><div class="detail-row"><span>${escape(c.confidence.name)}</span><b>${escape(c.confidence.value ?? "Not evaluated")}</b></div><small>${escape(c.confidence.scale)} · ${escape(c.confidence.rule)}</small><hr><div class="eyebrow">Missed core residues</div><div>${c.missed.length ? c.missed.map(key => residueButton(key, "missed")).join("") : "<small>None reported</small>"}</div><p class="compare-note">${escape(c.approval)}</p>`;
}

function residueButton(key, variant = "") {
  const label = manifest.signature.residues.find(r => r.key === key)?.label ?? key;
  return `<button class="residue-chip ${variant} ${state.residue === key ? "active" : ""}" data-residue="${escape(key)}" aria-pressed="${state.residue === key}">${escape(label)}</button>`;
}

function evidence() {
  return `<div class="section-head"><h2>Evidence before prediction.</h2><span class="tag">${manifest.evidence.length} evidence records</span></div><div class="cards"><article class="card"><div class="eyebrow">Disease hypothesis</div><h3>${escape(manifest.disease.name)}</h3><p>${escape(manifest.target.rationale)}</p><small>${escape(manifest.disease.id)}</small></article>${manifest.evidence.map(e => `<article class="card"><div class="eyebrow">${escape(e.id)}</div><h3>${escape(e.title)}</h3><p>${escape(e.claim)}</p>${e.url?.startsWith("https://") ? `<a href="${escape(e.url)}" target="_blank" rel="noopener noreferrer">Inspect source</a>` : "<small>No external citation supplied</small>"}</article>`).join("")}</div>`;
}

function siteDetails() {
  const structures = manifest.structures;
  const sig = manifest.signature;
  return `<div class="section-head"><h2>Structure is context, not a conclusion.</h2></div><div class="cards"><article class="card"><div class="eyebrow">Site assessment</div><h3>${escape(title(sig.addressability.status))}</h3><p>${escape(sig.addressability.reason)}</p></article>${structures.map(a => `<article class="card"><div class="eyebrow">${escape(a.origin)}</div><h3>${escape(a.id)}</h3><p>${escape(a.source)}</p><small>License: ${escape(a.license)}</small><p class="compare-note">${escape(bundle.display[a.id]?.alignment.reason ?? "No display alignment")}</p></article>`).join("")}</div>${bundle.warnings.map(w => `<div class="notice">${escape(w)}</div>`).join("")}`;
}

function signatureDetails() {
  const sig = manifest.signature;
  return `<div class="section-head"><h2>A signature built from the target side.</h2><span class="tag">${escape(shortSource[sig.source])}</span></div><div class="cards"><article class="card"><div class="eyebrow">Ensemble provenance</div><div class="stat">${sig.source === "boltzgen_consensus" ? `${sig.n_surviving} / ${sig.n_generated}` : "Reference"}</div><p>${sig.source === "boltzgen_consensus" ? "surviving / generated designs" : "Not a generated-design consensus"}</p><small>${escape(sig.reason)}</small></article><article class="card"><div class="eyebrow">Reading the surface</div><h3>Darker means more frequent engagement.</h3><p>Frequency has a fixed 0–1 scale. Core membership uses the upstream threshold; the viewer never retunes it.</p><small>${manifest.design_reference ? `Representative design: ${escape(manifest.design_reference.name)}. Protocol: ${escape(manifest.design_reference.protocol)}. Source: ${escape(manifest.design_reference.source)}. One pose is not the ensemble consensus.` : "No representative design was supplied. No reference binder is fabricated."}</small></article></div><div class="panel section table-scroll"><table><thead><tr><th>Residue</th><th>Engagement frequency</th><th>Core</th><th>Mapping</th></tr></thead><tbody>${sig.residues.map(r => `<tr><td>${residueButton(r.key)}</td><td>${Math.round(r.frequency * 100)}%<span class="bar"><span style="width:${r.frequency * 100}%"></span></span></td><td>${r.core ? "Core" : "Peripheral"}</td><td>${bundle.display[manifest.target.structure_id]?.residue_atoms[r.key] ? "Mapped" : "Unmapped — no highlight"}</td></tr>`).join("")}</tbody></table></div>`;
}

function filters() {
  return `<div class="filters"><label>Modality <select id="modality">${["small_molecule", "peptide", "biologic"].map(value => `<option value="${value}" ${state.modality === value ? "selected" : ""}>${title(value)}</option>`).join("")}</select></label><label>Novelty <select id="novelty">${[["all", "All pairings"], ["novel", "Novel only"], ["known", "Known only"]].map(([value, text]) => `<option value="${value}" ${state.novelty === value ? "selected" : ""}>${text}</option>`).join("")}</select></label><label><input type="checkbox" id="rejected" ${state.showRejected ? "checked" : ""}> Include rejected / unassessed</label></div>`;
}

function board() {
  const visible = rows(manifest, state);
  const isSuppressed = suppressed(manifest, state.modality);
  return `<div class="section-head"><h2>${state.chapter === 3 ? "Inspect an interface match." : "Ranked within modality."}</h2>${filters()}</div>${isSuppressed ? `<div class="notice">Small-molecule leaderboard suppressed. ${escape(manifest.signature.addressability.reason)} Records are inspection-only.</div>` : ""}<div class="panel table-scroll"><table><thead><tr><th>Rank</th><th>Molecule / pairing</th><th>Core coverage</th><th>Confidence gate</th><th>Decoy percentile</th></tr></thead><tbody>${visible.map(c => `<tr class="${c.id === state.candidate ? "selected" : ""}"><td>${!isSuppressed && c.rank !== null ? String(c.rank).padStart(2, "0") : "—"}</td><td><button class="table-button" data-candidate="${escape(c.id)}" aria-pressed="${c.id === state.candidate}">${escape(c.name)}<span class="subtext">${escape(title(c.novelty))}</span></button></td><td>${escape(percentage(c.coverage))}${c.coverage.value !== null ? `<span class="bar"><span style="width:${c.coverage.value * 100}%"></span></span>` : ""}</td><td><span class="tag ${c.confidence.gate === "passed" ? "" : "warning"}">${escape(title(c.confidence.gate))}</span></td><td>${c.decoys.status === "available" ? `${escape(c.decoys.percentile)}<span class="subtext">n = ${c.decoys.n}; within modality</span>` : "<span class='muted'>Uncalibrated</span>"}</td></tr>`).join("")}</tbody></table>${!visible.length ? "<div class='empty'>No candidates match these filters. No results have been filled in.</div>" : ""}</div><p class="subtle muted">Showing ${visible.length} of ${manifest.candidates.filter(c => c.modality === state.modality).length} records in this modality. Ordering supplied by the scientific pipeline: ${escape(visible[0]?.ranking_metric ?? "not available")}. Novel pairing means absent from the cited annotations, not proven novelty.</p>${heatmap(visible)}${candidate() ? `<div class="notice">${candidate().caveats.map(escape).join(" ") || "No candidate-specific caveats supplied."}</div>` : ""}`;
}

function heatmap(visible) {
  const core = manifest.signature.residues.filter(r => r.core);
  if (!core.length) return "<div class='notice'>No evaluable core signature. Coverage cannot be interpreted.</div>";
  return `<div class="section-head section"><h2>What is covered. What is missing.</h2><small>Select a residue to locate it in 3D.</small></div><div class="panel table-scroll"><table class="heatmap"><thead><tr><th>Candidate</th>${core.map(r => `<th><button data-residue="${escape(r.key)}">${escape(r.label)}<span class="subtext">${Math.round(r.frequency * 100)}% frequency</span></button></th>`).join("")}</tr></thead><tbody>${visible.map(c => `<tr><td><button class="table-button" data-candidate="${escape(c.id)}">${escape(c.name)}</button></td>${core.map(r => {const kind = c.engaged.includes(r.key) ? "engaged" : c.missed.includes(r.key) ? "missed" : "unavailable"; return `<td><button class="heat-cell ${kind} ${state.residue === r.key ? "active" : ""}" data-residue="${escape(r.key)}" aria-label="${escape(c.name)} / ${escape(r.label)}: ${kind}">${kind === "engaged" ? "+" : kind === "missed" ? "−" : "?"}</button></td>`;}).join("")}</tr>`).join("")}</tbody></table></div>`;
}

function nullChart(c) {
  if (!c || c.decoys.status !== "available") return `<div class="card"><div class="eyebrow">Null calibration</div><h3>Not evaluated</h3><p>${escape(c?.decoys.reason ?? "No candidate selected")}</p><small>No percentile or probability is invented.</small></div>`;
  const sorted = [...c.decoys.scores].sort((a, b) => a - b);
  const points = sorted.map((value, i) => `${40 + value * 500},${180 - ((i + 1) / sorted.length) * 140}`).join(" ");
  const x = c.coverage.value === null ? null : 40 + c.coverage.value * 500;
  return `<div class="card"><div class="eyebrow">Within-modality null / n = ${c.decoys.n}</div><h3>${escape(c.name)} · percentile ${escape(c.decoys.percentile)}</h3><svg class="chart" viewBox="0 0 600 220" role="img" aria-label="Decoy empirical cumulative distribution; candidate percentile ${escape(c.decoys.percentile)}"><path d="M40 30V180H555" fill="none" stroke="#94a3b8"/><polyline points="${points}" fill="none" stroke="#007f73" stroke-width="2"/>${x !== null ? `<path d="M${x} 35V180" stroke="#2563eb" stroke-width="2" stroke-dasharray="4 4"/>` : ""}<text x="40" y="201" font-size="11">0</text><text x="525" y="201" font-size="11">1</text><text x="200" y="218" font-size="11">Core coverage</text></svg><small>Empirical distribution, not a binding probability or significance test. Tie rule: ${escape(c.decoys.tie_rule)}</small></div>`;
}

function validation() {
  return `<div class="section-head"><h2>Keep the uncertainty in the picture.</h2></div><div class="cards">${nullChart(candidate())}${manifest.validation.map(v => `<article class="card"><div class="eyebrow">${escape(title(v.status))} / n = ${v.n}</div><h3>${escape(v.label)}</h3>${v.value !== null ? `<div class="stat">${escape(v.value)}</div>` : ""}<p>${escape(v.reason)}</p><small>Source: ${escape(v.source)}</small></article>`).join("")}</div><div class="notice">No superiority or significance claim is made by this viewer. Comparisons require upstream statistical analysis. Raw overlap without null calibration is explicitly uncalibrated.</div><section class="card section"><h3>Provenance</h3><dl class="provenance">${Object.entries(manifest.provenance).map(([key, value]) => `<dt>${escape(title(key))}</dt><dd>${escape(value)}</dd>`).join("")}</dl>${bundle.warnings.map(w => `<p>${escape(w)}</p>`).join("")}</section>`;
}

function render() {
  reconcile(manifest, state);
  const chapter = chapters[state.chapter];
  $("#hero").innerHTML = `<div class="hero"><div><div class="eyebrow">${escape(manifest.disease.name)} / ${chapter[1]}</div><h1>${chapter[2]}</h1><p>${chapter[3]}</p></div><div class="step">0${state.chapter + 1} <span class="muted">/ 06</span></div></div>`;
  document.querySelectorAll("[data-chapter]").forEach(el => {const active = Number(el.dataset.chapter) === state.chapter; el.classList.toggle("active", active); el.setAttribute("aria-current", active ? "step" : "false");});
  document.querySelectorAll("[data-view]").forEach(el => {const active = el.dataset.view === state.view; el.classList.toggle("active", active); el.setAttribute("aria-pressed", String(active));});
  $("#reference-label").textContent = state.chapter >= 2 && manifest.design_reference ? `${manifest.design_reference.name} (representative design)` : shortSource[manifest.signature.source];
  $("#candidate-label").textContent = candidate()?.name ?? "No candidate pose";
  $("#inspector").innerHTML = inspector();
  const frequency = state.view === "frequency" || state.chapter === 2;
  $("#legend").innerHTML = frequency ? '<span>Engagement frequency: light 0 → dark 1</span><span>Core membership uses the upstream threshold</span>' : `<span><i class="dot" style="background:#cbd5e1"></i>Target context</span><span><i class="dot" style="background:#007f73"></i>Covered core</span><span><i class="dot" style="background:#b45309"></i>Missed core</span><span><i class="dot" style="background:#2563eb"></i>Selected molecule</span>`;
  $("#chapter-content").innerHTML = [evidence, siteDetails, signatureDetails, board, board, validation][state.chapter]();
  $("#previous").disabled = state.chapter === 0;
  $("#next").disabled = state.chapter === 5;
  $("#residue-selection").textContent = state.residue ? `Selected residue: ${manifest.signature.residues.find(r => r.key === state.residue)?.label ?? state.residue}` : "";
  $("#modality")?.addEventListener("change", event => {state.modality = event.target.value; render();});
  $("#novelty")?.addEventListener("change", event => {state.novelty = event.target.value; render();});
  $("#rejected")?.addEventListener("change", event => {state.showRejected = event.target.checked; render();});
  stage.draw(bundle, state).catch(() => {$("#viewer-status").textContent = "Structure rendering failed; tabular explanations remain available.";});
}

try {
  const response = await fetch("bundle.json");
  if (!response.ok) throw new Error(`Bundle request failed (${response.status})`);
  bundle = await response.json();
  manifest = bundle.manifest;
  if (manifest?.schema_version !== "1.0.0" || !Array.isArray(manifest.candidates) || !bundle.display) throw new Error("Unsupported or invalid bundle; rebuild with the visualization validator.");
  state = initialState(manifest);
  shell();
  render();
} catch (error) {
  $("#app").innerHTML = `<main class="loading"><h1>Snapshot unavailable</h1><p>${escape(error.message)}</p><p>Build a validated export and serve that directory over localhost. The interactive viewer does not fetch scientific data from the network.</p></main>`;
}
