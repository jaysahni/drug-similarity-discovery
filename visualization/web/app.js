import {initialState, percentage, reconcile, rows, suppressed} from "./state.js";
import {MolecularStage} from "./viewer.js";
import {renderPipeline} from "./pipeline.js";

const $ = selector => document.querySelector(selector);
const escape = value => String(value ?? "").replace(/[&<>"']/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));
const title = value => value.replaceAll("_", " ");
const shortSource = {boltzgen_consensus: "BoltzGen consensus", p2rank_geometry: "Pocket geometry", known_ligand: "Known-ligand reference", hybrid: "Hybrid signature"};
let bundle, manifest, state, stage, activeTab = "candidates", opened = false;
const candidate = () => manifest.candidates.find(c => c.id === state.candidate);
const isDemo = () => manifest.data_origin === "synthetic_fixture";
const isReference = () => manifest.data_origin === "reference_example";

function shell() {
  $("#app").innerHTML = `
    <aside class="sidebar">
      <a class="brand" href="./"><span class="brand-mark">a</span> autorepurpose<span class="brand-period">.</span></a>
      <button class="new-chat" id="new-chat"><span>＋</span> New chat</button>
      <div class="sidebar-section">Workspace</div>
      <button class="workspace-item active" id="explore"><span>◌</span> Research</button>
      <div class="sidebar-section saved-label">Saved results</div>
      <button class="snapshot" data-open><span class="snapshot-dot"></span><span>${escape(isDemo() ? "Example results" : manifest.disease.name)}<small>${escape(isDemo() ? "Synthetic demo" : manifest.target.name)}</small></span></button>
      <div class="sidebar-bottom"><span class="avatar">R</span><div>Research workspace<small>Local preview</small></div><span class="online-dot"></span></div>
    </aside>
    <div class="main-wrap"><header class="topbar"><span>Research <span class="slash">/</span> <span class="muted" id="thread-title">New chat</span></span><span class="preview-pill"><span></span> Preview</span></header>
    <main id="main" class="content" tabindex="-1">
      <section class="welcome" id="welcome"><h1>Disease Research</h1><p>Enter a disease to get started.</p></section>
      <section id="conversation" class="conversation" aria-label="Research conversation" aria-live="polite"></section>
      <section id="result" class="result" hidden>
        <div class="assistant-heading"><span class="assistant-avatar">a</span><b>AutoRepurpose</b><span id="origin-banner" class="tag"></span></div>
        <div class="response-intro"><h2 id="result-title"></h2><p id="result-description"></p></div>
        <div id="pipeline-root"></div>
        <div class="result-toolbar"><div class="tabs" role="tablist" aria-label="Research results">${[["candidates","Candidates"],["structure","Structure"],["sources","Sources"]].map(([id,label])=>`<button role="tab" id="tab-${id}" aria-controls="result-panel" data-tab="${id}" aria-selected="${id === 'candidates'}">${label}</button>`).join("")}</div><a href="report.html" target="_blank" rel="noopener" class="report-link">Full report</a></div>
        <div id="result-panel" role="tabpanel" aria-labelledby="tab-candidates"><div id="chapter-content"></div>
        <div id="molecular-section" hidden><div class="panel"><div class="panel-heading"><h2>Molecular structure</h2><div class="segmented" role="group" aria-label="Molecular view">${[["candidate","Candidate"],["binding_site","Candidate with binding site"]].map(([id,label])=>`<button data-view="${id}">${label}</button>`).join("")}</div></div><div class="molecular-label"><b id="candidate-label"></b></div><div id="viewer-candidate" class="viewer" role="img" aria-label="Selected candidate structure"></div><div id="legend" class="legend"></div><div id="viewer-status" class="viewer-status" role="status"></div></div><div id="residue-selection" class="subtle" role="status"></div></div>
        <aside id="inspector" class="inspector" aria-label="Selected candidate"></aside></div>
        <details class="detail-disclosure" id="validation"><summary>Validation & run details</summary><div id="validation-content"></div><a href="bundle.json" download>Download results</a></details>
      </section>
      <div class="composer-wrap" id="composer-wrap"><form class="composer" id="chat-form"><label class="screen-reader" for="prompt">Disease or research question</label><textarea id="prompt" rows="2" maxlength="2000" placeholder="Enter a disease or question" required></textarea><div class="composer-bottom"><span class="composer-mode">Research</span><button type="submit" class="send" aria-label="Send research question">↑</button></div></form><p class="composer-note">Preview only. Live search and analysis are not connected.</p></div>
      <section class="starter-content" id="starter-content"><div class="suggestions"><button data-prompt="Colorectal cancer">Colorectal cancer</button><button data-prompt="Alzheimer’s disease">Alzheimer’s disease</button><button data-prompt="Type 2 diabetes">Type 2 diabetes</button></div><div class="workflow"><div><span>01</span><b>Targets</b><p>Literature review</p></div><div><span>02</span><b>Binders</b><p>Molecular design</p></div><div><span>03</span><b>Drug comparison</b><p>Similarity and coverage</p></div><div><span>04</span><b>Candidates</b><p>Results for review</p></div></div><button class="example-link" data-open>View ${isDemo() ? "example results" : "saved results"}</button></section>

    </main></div>`;
  $("#chat-form").addEventListener("submit", submitPrompt);
  $("#prompt").addEventListener("keydown", event => {if (event.key === "Enter" && !event.shiftKey && !event.isComposing) {event.preventDefault(); $("#chat-form").requestSubmit();}});
  $("#new-chat").addEventListener("click", reset);
  $("#explore").addEventListener("click", reset);
  document.addEventListener("click", event => {
    const button = event.target.closest("button");
    if (!button) return;
    if (button.hasAttribute("data-open")) openSnapshot();
    if (button.dataset.prompt) {$("#prompt").value = button.dataset.prompt; $("#prompt").focus();}
    if (button.dataset.tab) {activeTab = button.dataset.tab; renderResults();}
    if (button.dataset.candidate) {state.candidate = button.dataset.candidate; state.view = "candidate"; renderResults();}
    if (button.dataset.residue) {state.view = "binding_site"; state.residue = state.residue === button.dataset.residue ? null : button.dataset.residue; activeTab = "structure"; renderResults();}
    if (button.dataset.view) {state.view = button.dataset.view; renderResults();}
  });
  document.querySelector('.tabs').addEventListener('keydown', event => {
    if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
    const tabs = [...document.querySelectorAll('[data-tab]')];
    const current = tabs.indexOf(document.activeElement);
    if (current < 0) return;
    event.preventDefault();
    const next = event.key === 'Home' ? 0 : event.key === 'End' ? tabs.length-1 : (current + (event.key === 'ArrowRight' ? 1 : -1) + tabs.length) % tabs.length;
    tabs[next].click(); tabs[next].focus();
  });
}

function begin() {
  $("#welcome").hidden = true;
  $("#starter-content").hidden = true;
  $("#main").classList.add("has-conversation");
}
function reset() {
  opened = false;
  $("#conversation").replaceChildren();
  $("#result").hidden = true;
  $("#welcome").hidden = false;
  $("#starter-content").hidden = false;
  $("#main").classList.remove("has-conversation");
  $("#thread-title").textContent = "New chat";
  $("#prompt").value = "";
  $("#prompt").focus();
}
function submitPrompt(event) {
  event.preventDefault();
  const prompt = $("#prompt").value.trim();
  if (!prompt) return;
  begin();
  $("#result").hidden = true;
  opened = false;
  const message = document.createElement("div");
  message.className = "exchange";
  message.innerHTML = `<div class="user-message">${escape(prompt)}</div><div class="assistant-heading"><span class="assistant-avatar">a</span><b>AutoRepurpose</b></div><p class="assistant-message">New disease analysis isn’t connected yet. ${isDemo() ? "Example results use synthetic data." : `A saved analysis for <strong>${escape(manifest.disease.name)}</strong> is available.`}</p><button class="inline-action" data-open>${isDemo() ? "Open example results" : "Open saved analysis"}</button>`;
  $("#conversation").append(message);
  $("#thread-title").textContent = prompt.length > 45 ? prompt.slice(0,45)+"…" : prompt;
  $("#prompt").value = "";
  message.scrollIntoView({block:"start",behavior:"smooth"});
}
function openSnapshot() {
  begin();
  if (!opened) {
    const message = document.createElement("div");
    message.className = "user-message";
    message.textContent = isDemo() ? "Show me the example results." : `Open the saved ${manifest.disease.name} analysis.`;
    $("#conversation").append(message);
  }
  opened = true;
  state = initialState(manifest);
  activeTab = isReference() ? "structure" : "candidates";
  $("#result").hidden = false;
  $("#thread-title").textContent = isDemo() ? "Example results" : manifest.disease.name;
  $("#origin-banner").textContent = isDemo() ? "Synthetic example" : manifest.data_origin === "computed" ? "Saved analysis" : "Reference example";
  $("#result-title").textContent = isDemo() ? "Example results" : `${manifest.disease.name} · ${manifest.target.name}`;
  $("#result-description").textContent = isDemo() ? "Synthetic molecules and scores for testing." : isReference() ? "" : `${manifest.candidates.length} candidate records in this saved analysis.`;
  $("#result-description").hidden = isReference();
  renderPipeline($("#pipeline-root"), manifest, (tab, view) => {
    activeTab = tab;
    if (view) state.view = view;
    renderResults();
    $("#result-panel").scrollIntoView({block: "nearest", behavior: "smooth"});
  });
  renderResults();
  $("#result").scrollIntoView({block:"start",behavior:"smooth"});
}
function renderResults() {
  reconcile(manifest,state);
  state.chapter = 3;
  document.querySelectorAll('[data-tab]').forEach(el=>{const active=el.dataset.tab===activeTab; el.setAttribute('aria-selected',String(active)); el.tabIndex=active ? 0 : -1;});
  $("#result-panel").setAttribute("aria-labelledby",`tab-${activeTab}`);
  $("#chapter-content").innerHTML = activeTab === "candidates" ? board() : activeTab === "sources" ? evidence() : "";
  $("#molecular-section").hidden = activeTab !== "structure";
  $("#inspector").hidden = activeTab === "sources";
  $("#inspector").innerHTML = inspector();
  $("#validation-content").innerHTML = validation();
  $("#modality")?.addEventListener("change", e=>{state.modality=e.target.value;renderResults();});
  $("#novelty")?.addEventListener("change", e=>{state.novelty=e.target.value;renderResults();});
  $("#rejected")?.addEventListener("change", e=>{state.showRejected=e.target.checked;renderResults();});
  if (activeTab === "structure") {
    if (!stage) stage = new MolecularStage($("#viewer-candidate"),$("#viewer-status"),key=>{state.residue=state.residue===key?null:key;renderResults();});
    $("#candidate-label").textContent = candidate()?.name ?? "No candidate pose";
    $("#legend").innerHTML = '<span><i class="dot blue"></i>Candidate</span>' + (state.view === "binding_site" ? '<span><i class="dot teal"></i>Reference binding site</span><span>Target protein in grey</span>' : "");
    $("#residue-selection").textContent = state.residue ? `Selected residue: ${manifest.signature.residues.find(r=>r.key===state.residue)?.label ?? state.residue}` : "";
    document.querySelectorAll('[data-view]').forEach(el=>{el.classList.toggle('active',el.dataset.view===state.view);el.setAttribute('aria-pressed',String(el.dataset.view===state.view));});
    stage.draw(bundle,state).catch(()=>{$("#viewer-status").textContent="Structure unavailable. Saved scores are still accessible.";});
  }
}

function inspector() {
  const c = candidate();
  const sig = manifest.signature;
  const core = sig.residues.filter(r => r.core);
  if (state.chapter < 3) {
    return `<div class="eyebrow">${state.chapter === 0 ? "The starting point" : "The interface"}</div><h2>${escape(manifest.target.name)}</h2><p>${escape(manifest.target.site_kind)} · ${escape(manifest.target.id)}</p><div class="score">${sig.status === "available" ? core.length : "—"}</div><div class="score-label">${sig.status === "available" ? "consensus-defined core residues" : "Signature not available"}</div><hr><div class="detail-row"><span>Signature source</span><b>${escape(shortSource[sig.source])}</b></div><div class="detail-row"><span>Core threshold</span><b>${Math.round(sig.threshold * 100)}%</b></div><div class="detail-row"><span>Addressability</span><b>${escape(title(sig.addressability.status))}</b></div><hr><p>${escape(sig.addressability.reason)}</p><small>Consensus contact frequency is not evidence of therapeutic efficacy.</small>`;
  }
  if (!c) return `<div class="eyebrow">Selected candidate</div><h2>No candidate selected</h2><p>No candidates match the current filters.</p>`;
  if (isReference()) return `<div class="eyebrow">Experimental reference</div><h2>${escape(c.name)}</h2><p>${escape(manifest.provenance.method ?? "Recorded structure")}</p><p>${escape(c.approval)}</p>${c.source.startsWith("https://") ? `<a href="${escape(c.source)}" target="_blank" rel="noopener">View structure source</a>` : ""}`;
  const coverage = percentage(c.coverage);
  return `<div class="eyebrow">Selected candidate</div><h2>${escape(c.name)}</h2><span class="tag ${c.confidence.gate === "passed" ? "" : "warning"}">${escape(title(c.novelty))}</span><div class="score">${escape(coverage)}</div><div class="score-label">${c.coverage.status === "available" ? `${c.engaged.length} of ${core.length} core residues engaged` : escape(c.coverage.reason)}</div><hr><div class="detail-row"><span>Confidence gate</span><b>${escape(title(c.confidence.gate))}</b></div><div class="detail-row"><span>${escape(c.confidence.name)}</span><b>${escape(c.confidence.value ?? "Not evaluated")}</b></div><small>${escape(c.confidence.scale)} · ${escape(c.confidence.rule)}</small><hr><div class="eyebrow">Missed core residues</div><div>${c.missed.length ? c.missed.map(key => residueButton(key, "missed")).join("") : "<small>None reported</small>"}</div><p class="compare-note">${escape(c.approval)}</p>`;
}

function residueButton(key, variant = "") {
  const label = manifest.signature.residues.find(r => r.key === key)?.label ?? key;
  return `<button class="residue-chip ${variant} ${state.residue === key ? "active" : ""}" data-residue="${escape(key)}" aria-pressed="${state.residue === key}">${escape(label)}</button>`;
}

function evidence() {
  return `<div class="section-head"><h2>Sources</h2><span class="tag">${manifest.evidence.length} evidence records</span></div><div class="cards"><article class="card"><div class="eyebrow">Disease</div><h3>${escape(manifest.disease.name)}</h3><p>${escape(manifest.target.rationale)}</p><small>${escape(manifest.disease.id)}</small></article>${manifest.evidence.map(e => `<article class="card"><div class="eyebrow">${escape(e.id)}</div><h3>${escape(e.title)}</h3><p>${escape(e.claim)}</p>${e.url?.startsWith("https://") ? `<a href="${escape(e.url)}" target="_blank" rel="noopener noreferrer">View source</a>` : "<small>No external citation supplied</small>"}</article>`).join("")}</div>`;
}

function siteDetails() {
  const structures = manifest.structures;
  const sig = manifest.signature;
  return `<div class="section-head"><h2>Structure details</h2></div><div class="cards"><article class="card"><div class="eyebrow">Site assessment</div><h3>${escape(title(sig.addressability.status))}</h3><p>${escape(sig.addressability.reason)}</p></article>${structures.map(a => `<article class="card"><div class="eyebrow">${escape(a.origin)}</div><h3>${escape(a.id)}</h3><p>${escape(a.source)}</p><small>License: ${escape(a.license)}</small><p class="compare-note">${escape(bundle.display[a.id]?.alignment.reason ?? "No display alignment")}</p></article>`).join("")}</div>${bundle.warnings.map(w => `<div class="notice">${escape(w)}</div>`).join("")}`;
}

function signatureDetails() {
  const sig = manifest.signature;
  return `<div class="section-head"><h2>Contact signature</h2><span class="tag">${escape(shortSource[sig.source])}</span></div><div class="cards"><article class="card"><div class="eyebrow">Ensemble provenance</div><div class="stat">${sig.source === "boltzgen_consensus" ? `${sig.n_surviving} / ${sig.n_generated}` : "Reference"}</div><p>${sig.source === "boltzgen_consensus" ? "surviving / generated designs" : "Not a generated-design consensus"}</p><small>${escape(sig.reason)}</small></article><article class="card"><div class="eyebrow">Contact frequency</div><h3>Darker means more frequent engagement.</h3><p>Frequency has a fixed 0–1 scale. Core membership uses the upstream threshold; the viewer never retunes it.</p><small>${manifest.design_reference ? `Representative design: ${escape(manifest.design_reference.name)}. Protocol: ${escape(manifest.design_reference.protocol)}. Source: ${escape(manifest.design_reference.source)}. One pose is not the ensemble consensus.` : "No representative design supplied."}</small></article></div><div class="panel section table-scroll"><table><thead><tr><th>Residue</th><th>Engagement frequency</th><th>Core</th><th>Mapping</th></tr></thead><tbody>${sig.residues.map(r => `<tr><td>${residueButton(r.key)}</td><td>${Math.round(r.frequency * 100)}%<span class="bar"><span style="width:${r.frequency * 100}%"></span></span></td><td>${r.core ? "Core" : "Peripheral"}</td><td>${bundle.display[manifest.target.structure_id]?.residue_atoms[r.key] ? "Mapped" : "Unmapped — no highlight"}</td></tr>`).join("")}</tbody></table></div>`;
}

function filters() {
  return `<div class="filters"><label>Modality <select id="modality">${["small_molecule", "peptide", "biologic"].map(value => `<option value="${value}" ${state.modality === value ? "selected" : ""}>${title(value)}</option>`).join("")}</select></label><label>Novelty <select id="novelty">${[["all", "All pairings"], ["novel", "Novel only"], ["known", "Known only"]].map(([value, text]) => `<option value="${value}" ${state.novelty === value ? "selected" : ""}>${text}</option>`).join("")}</select></label><label><input type="checkbox" id="rejected" ${state.showRejected ? "checked" : ""}> Include rejected / unassessed</label></div>`;
}

function board() {
  const visible = rows(manifest, state);
  const isSuppressed = suppressed(manifest, state.modality);
  return `<div class="section-head"><h2>Candidates</h2>${filters()}</div>${isSuppressed ? `<div class="notice">Small-molecule ranking unavailable. ${escape(manifest.signature.addressability.reason)}</div>` : ""}<div class="panel table-scroll"><table><thead><tr><th>Rank</th><th>Molecule / pairing</th><th>Core coverage</th><th>Confidence gate</th><th>Decoy percentile</th></tr></thead><tbody>${visible.map(c => `<tr class="${c.id === state.candidate ? "selected" : ""}"><td>${!isSuppressed && c.rank !== null ? String(c.rank).padStart(2, "0") : "—"}</td><td><button class="table-button" data-candidate="${escape(c.id)}" aria-pressed="${c.id === state.candidate}">${escape(c.name)}<span class="subtext">${escape(title(c.novelty))}</span></button></td><td>${escape(percentage(c.coverage))}${c.coverage.value !== null ? `<span class="bar"><span style="width:${c.coverage.value * 100}%"></span></span>` : ""}</td><td><span class="tag ${c.confidence.gate === "passed" ? "" : "warning"}">${escape(title(c.confidence.gate))}</span></td><td>${c.decoys.status === "available" ? `${escape(c.decoys.percentile)}<span class="subtext">n = ${c.decoys.n}; within modality</span>` : "<span class='muted'>Uncalibrated</span>"}</td></tr>`).join("")}</tbody></table>${!visible.length ? "<div class='empty'>No candidates match these filters.</div>" : ""}</div><p class="subtle muted">Showing ${visible.length} of ${manifest.candidates.filter(c => c.modality === state.modality).length} records in this modality. Ordering supplied by the scientific pipeline: ${escape(visible[0]?.ranking_metric ?? "not available")}. Novel pairing means absent from the cited annotations, not proven novelty.</p><details class="detail-disclosure"><summary>View residue coverage</summary>${heatmap(visible)}</details>${candidate() ? `<details class="detail-disclosure"><summary>Candidate notes</summary><div class="notice">${candidate().caveats.map(escape).join(" ") || "No candidate-specific caveats supplied."}</div></details>` : ""}`;
}

function heatmap(visible) {
  const core = manifest.signature.residues.filter(r => r.core);
  if (!core.length) return "<div class='notice'>No evaluable core signature. Coverage cannot be interpreted.</div>";
  return `<div class="section-head section"><h2>Residue coverage</h2><small>Select a residue to locate it in 3D.</small></div><div class="panel table-scroll"><table class="heatmap"><thead><tr><th>Candidate</th>${core.map(r => `<th><button data-residue="${escape(r.key)}">${escape(r.label)}<span class="subtext">${Math.round(r.frequency * 100)}% frequency</span></button></th>`).join("")}</tr></thead><tbody>${visible.map(c => `<tr><td><button class="table-button" data-candidate="${escape(c.id)}">${escape(c.name)}</button></td>${core.map(r => {const kind = c.engaged.includes(r.key) ? "engaged" : c.missed.includes(r.key) ? "missed" : "unavailable"; return `<td><button class="heat-cell ${kind} ${state.residue === r.key ? "active" : ""}" data-residue="${escape(r.key)}" aria-label="${escape(c.name)} / ${escape(r.label)}: ${kind}">${kind === "engaged" ? "+" : kind === "missed" ? "−" : "?"}</button></td>`;}).join("")}</tr>`).join("")}</tbody></table></div>`;
}

function nullChart(c) {
  if (!c || c.decoys.status !== "available") return `<article class="card"><h3>Calibration unavailable</h3><p>${escape(c?.decoys.reason ?? "No candidate selected")}</p></article>`;
  return `<article class="card"><h3>Decoy percentile · ${escape(c.decoys.percentile)}</h3><meter min="0" max="100" value="${c.decoys.percentile}">${c.decoys.percentile}</meter><p>${c.decoys.n} within-modality comparisons · ${escape(c.decoys.tie_rule)}</p><small>Not a binding probability.</small></article>`;
}

function validation() {
  return `<div class="section-head"><h2>Validation</h2></div><div class="cards">${nullChart(candidate())}${manifest.validation.map(v => `<article class="card"><div class="eyebrow">${escape(title(v.status))} / n = ${v.n}</div><h3>${escape(v.label)}</h3>${v.value !== null ? `<div class="stat">${escape(v.value)}</div>` : ""}<p>${escape(v.reason)}</p><small>Source: ${escape(v.source)}</small></article>`).join("")}</div><div class="notice">Coverage scores alone do not establish statistical significance.</div><section class="card section"><h3>Data sources</h3><dl class="provenance">${Object.entries(manifest.provenance).map(([key, value]) => `<dt>${escape(title(key))}</dt><dd>${escape(value)}</dd>`).join("")}</dl>${bundle.warnings.map(w => `<p>${escape(w)}</p>`).join("")}</section>`;
}


try {
  const response = await fetch("bundle.json");
  if (!response.ok) throw new Error("Could not load the saved analysis.");
  bundle = await response.json();
  manifest = bundle.manifest;
  if (manifest?.schema_version !== "1.0.0" || !Array.isArray(manifest.candidates) || !bundle.display) throw new Error("This saved analysis needs to be rebuilt.");
  state = initialState(manifest);
  shell();
} catch (error) {
  $("#app").innerHTML = `<main class="loading"><h1>Preview unavailable</h1><p>${escape(error.message)}</p><p>Build a visualization export and open it through the local server.</p></main>`;
}
