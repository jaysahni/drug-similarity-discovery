import {initialState, percentage, reconcile, rows, suppressed} from "./state.js";
import {MolecularStage} from "./viewer.js";
import {renderPipeline} from "./pipeline.js";
import {newChat, addFolder, updateChat, removeChat, removeFolder, findSavedResult, saveResult, removeSavedResult, readWorkspace, writeWorkspace} from "./chats.js";

const $ = selector => document.querySelector(selector);
const escape = value => String(value ?? "").replace(/[&<>"']/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));
const title = value => value.replaceAll("_", " ");
const shortSource = {boltzgen_consensus: "BoltzGen consensus", p2rank_geometry: "Pocket geometry", known_ligand: "Known-ligand reference", hybrid: "Hybrid signature"};
let bundle, manifest, state, stage, activeStep = "candidates", opened = false;
let workspace, storage, storageKey, storageWritable = true, organizerTarget = null;
const currentChat = () => workspace.chats.find(c => c.id === workspace.activeId);
const candidate = () => manifest.candidates.find(c => c.id === state.candidate);
const isDemo = () => manifest.data_origin === "synthetic_fixture";
const isReference = () => manifest.data_origin === "reference_example";

function shell() {
  $("#app").innerHTML = `
    <aside class="sidebar">
      <a class="brand" href="./"><span class="brand-mark">r</span> Rebind</a>
      <button class="new-chat" id="new-chat"><span>＋</span> New chat</button>
      <div class="chat-section-heading"><span>Chats</span><button id="new-folder" aria-label="New folder" title="New folder">＋ Folder</button></div>
      <div id="chat-list" class="chat-list" role="tablist" aria-label="Chats" aria-orientation="vertical"></div>
      <section class="saved-results-section" aria-labelledby="saved-results-heading"><h2 id="saved-results-heading">Saved results</h2><div id="saved-results-list"></div></section>
      <div class="sidebar-bottom"><span class="avatar">R</span><div>Research workspace<small id="storage-status" role="status">Saved in this browser</small></div><span class="online-dot"></span></div>
    </aside>
    <div class="main-wrap"><header class="topbar"><span>Research <span class="slash">/</span> <span class="muted" id="thread-title">New chat</span></span><span class="preview-pill"><span></span> Preview</span></header>
    <dialog id="organizer" aria-labelledby="organizer-title"><form id="organizer-form"><h2 id="organizer-title"></h2><label for="organizer-name">Name</label><input id="organizer-name" maxlength="80" required autocomplete="off"><label id="folder-field" for="organizer-folder">Folder<select id="organizer-folder"></select></label><p id="organizer-note"></p><p id="organizer-error" role="alert"></p><div class="organizer-actions"><button type="button" id="organizer-delete" class="delete-action">Delete</button><button type="button" id="organizer-cancel">Cancel</button><button type="submit" class="primary">Save</button></div></form></dialog>
    <main id="main" class="content" tabindex="-1"><div id="chat-panel" role="tabpanel">
      <section class="welcome" id="welcome"><h1>Disease Research</h1><p>Enter a disease to get started.</p></section>
      <section id="conversation" class="conversation" aria-label="Research conversation" aria-live="polite"></section>
      <section id="result" class="result" hidden>
        <div class="assistant-heading"><span class="assistant-avatar">r</span><b>Rebind</b><span id="origin-banner" class="tag"></span></div>
        <div class="response-intro"><div class="result-heading"><h2 id="result-title"></h2><button id="save-result" type="button" aria-pressed="false">Save result</button></div><p id="result-description"></p></div>
        <div id="pipeline-root"></div>
        <div id="result-panel" role="tabpanel" aria-labelledby="step-candidates"><div id="chapter-content"></div>
        <div id="molecular-section" hidden><div class="panel"><div class="panel-heading"><h2>Molecular structure</h2><div class="segmented" role="group" aria-label="Molecular view">${[["candidate","Candidate"],["binding_site","Candidate with binding site"]].map(([id,label])=>`<button data-view="${id}">${label}</button>`).join("")}</div></div><label class="molecular-label" for="structure-candidate">Molecule <select id="structure-candidate"></select></label><div id="viewer-candidate" class="viewer" role="img" aria-label="Selected candidate structure"></div><div id="legend" class="legend"></div><div id="viewer-status" class="viewer-status" role="status"></div></div><div id="residue-selection" class="subtle" role="status"></div></div>
        </div>
        <details class="detail-disclosure" id="validation"><summary>Validation & run details</summary><div id="validation-content"></div></details><div class="result-actions"><a href="report.html" target="_blank" rel="noopener">Full report</a><a href="bundle.json" download>Download results</a></div>
      </section>
      <div class="composer-wrap" id="composer-wrap"><form class="composer" id="chat-form"><label class="screen-reader" for="prompt">Disease or research question</label><textarea id="prompt" rows="2" maxlength="2000" placeholder="Enter a disease or question" required></textarea><div class="composer-bottom"><button type="submit" class="send" aria-label="Send research question">↑</button></div></form><p class="composer-note">Preview only. Live search and analysis are not connected.</p></div>


    </div></main></div>`;
  $("#chat-form").addEventListener("submit", submitPrompt);
  $("#prompt").addEventListener("keydown", event => {if (event.key === "Enter" && !event.shiftKey && !event.isComposing) {event.preventDefault(); $("#chat-form").requestSubmit();}});
  $("#new-chat").addEventListener("click", () => startChat());
  $("#save-result").addEventListener("click", toggleSavedResult);
  $("#new-folder").addEventListener("click", () => openOrganizer('new-folder'));
  $("#organizer-form").addEventListener("submit", saveOrganizer);
  $("#organizer-cancel").onclick = () => $("#organizer").close();
  $("#organizer-delete").onclick = deleteOrganizedItem;
  $("#prompt").addEventListener("input", () => {currentChat().draft = $("#prompt").value; persistChats();});
  $("#chat-list").addEventListener("keydown", event => {
    if (!event.target.matches('[data-chat]') || !['ArrowUp','ArrowDown','Home','End'].includes(event.key)) return;
    event.preventDefault();
    const tabs = [...document.querySelectorAll('[data-chat]')].filter(tab => tab.getClientRects().length);
    const index = tabs.indexOf(event.target);
    const next = event.key === 'Home' ? 0 : event.key === 'End' ? tabs.length-1 : (index + (event.key === 'ArrowUp' ? -1 : 1) + tabs.length) % tabs.length;
    const chatId = tabs[next].dataset.chat;
    switchChat(chatId);
    document.querySelector(`[data-chat="${chatId}"]`).focus();
  });
  document.addEventListener("click", event => {
    const button = event.target.closest("button");
    if (!button) return;
    if (button.dataset.savedResult) openSavedResult(button.dataset.savedResult);
    if (button.dataset.removeSaved) {
      removeSavedResult(workspace, button.dataset.removeSaved);
      persistChats(); renderSavedResults();
    }
    if (button.dataset.chat) switchChat(button.dataset.chat);
    if (button.hasAttribute('data-new-in-folder')) startChat(button.dataset.newInFolder);
    if (button.dataset.editChat) openOrganizer('chat', button.dataset.editChat);
    if (button.dataset.editFolder) openOrganizer('folder', button.dataset.editFolder);
    if (button.dataset.toggleFolder) {
      const folder = workspace.folders.find(f => f.id === button.dataset.toggleFolder);
      folder.collapsed = !folder.collapsed; persistChats(); renderChatList();
    }
    if (button.hasAttribute("data-open")) openSnapshot();
    if (button.dataset.candidate) {state.candidate = button.dataset.candidate; state.view = "candidate"; renderResults();}
    if (button.dataset.residue) {state.view = "binding_site"; state.residue = state.residue === button.dataset.residue ? null : button.dataset.residue; activeStep = "signature"; renderResults();}
    if (button.dataset.view) {state.view = button.dataset.view; renderResults();}
  });

}

function begin() {
  $("#welcome").hidden = true;
  $("#main").classList.add("has-conversation");
}
function persistChats() {
  const saved = storageWritable && writeWorkspace(storage, storageKey, workspace);
  if ($("#storage-status")) $("#storage-status").textContent = saved ? "Saved in this browser" : "Not saved · browser storage unavailable";
}
function rememberChat() {
  const chat = currentChat();
  chat.draft = $("#prompt").value;
  chat.snapshotOpen = opened;
  chat.activeStep = activeStep;
  chat.viewState = {...state};
  persistChats();
}
function renderChatList() {
  const row = chat => `<div class="chat-row ${chat.id === workspace.activeId ? 'active' : ''}"><button role="tab" id="chat-${chat.id}" data-chat="${chat.id}" aria-controls="chat-panel" aria-selected="${chat.id === workspace.activeId}" tabindex="${chat.id === workspace.activeId ? 0 : -1}" title="${escape(chat.name)}"><span>${escape(chat.name)}</span></button><button class="chat-menu" data-edit-chat="${chat.id}" aria-label="Manage ${escape(chat.name)}" title="Rename or move chat">⋯</button></div>`;
  $("#chat-list").innerHTML = workspace.folders.map(folder => {
    const chats = workspace.chats.filter(c => c.folderId === folder.id);
    return `<section class="chat-folder"><div class="folder-heading"><button data-toggle-folder="${folder.id}" aria-expanded="${!folder.collapsed}" aria-controls="folder-${folder.id}"><span aria-hidden="true">${folder.collapsed ? '▸' : '▾'}</span><span class="folder-name">${escape(folder.name)}</span><small>${chats.length}</small></button><button data-new-in-folder="${folder.id}" aria-label="New chat in ${escape(folder.name)}" title="New chat in folder">＋</button><button data-edit-folder="${folder.id}" aria-label="Manage folder ${escape(folder.name)}">⋯</button></div><div id="folder-${folder.id}" ${folder.collapsed ? 'hidden' : ''}>${chats.map(row).join('') || '<p class="folder-empty">No chats yet</p>'}</div></section>`;
  }).join('') + `<div class="unfiled-chats">${workspace.folders.length ? '<div class="unfiled-label">Unfiled</div>' : ''}${workspace.chats.filter(c => !c.folderId).map(row).join('')}</div>`;
  $("#chat-panel").setAttribute('aria-labelledby', `chat-${workspace.activeId}`);
  renderSavedResults();
}
function renderSavedResults() {
  $("#saved-results-list").innerHTML = workspace.savedResults.map(result => `<div class="saved-result-row"><button class="saved-result-open" data-saved-result="${result.id}" title="${escape(result.name)}"><span>${escape(result.name)}</span><small>${escape(result.subtitle)}</small></button><button class="saved-result-remove" data-remove-saved="${result.id}" aria-label="Remove ${escape(result.name)} from saved results" title="Remove from saved results">×</button></div>`).join('') || '<p class="saved-results-empty">No saved results yet.</p>';
  const saved = findSavedResult(workspace, manifest.run_id, candidate()?.id ?? null);
  $("#save-result").textContent = saved ? 'Saved' : 'Save result';
  $("#save-result").setAttribute('aria-pressed', String(Boolean(saved)));
  $("#save-result").title = saved ? 'Remove from saved results' : 'Save this result and its current view';
}
function toggleSavedResult() {
  if (!opened) return;
  const c = candidate();
  const saved = findSavedResult(workspace, manifest.run_id, c?.id ?? null);
  if (saved) removeSavedResult(workspace, saved.id);
  else saveResult(workspace, {
    runId: manifest.run_id, candidateId: c?.id ?? null,
    name: isDemo() ? 'Example results' : manifest.disease.name,
    subtitle: c?.name ?? manifest.target.name, chatId: workspace.activeId,
    viewState: state, activeStep,
  });
  persistChats(); renderSavedResults();
}
function openSavedResult(resultId) {
  const result = workspace.savedResults.find(item => item.id === resultId);
  if (!result || result.runId !== manifest.run_id) return;
  rememberChat();
  let chat = workspace.chats.find(item => item.id === result.chatId);
  if (!chat) {
    chat = newChat(workspace);
    chat.name = result.name; chat.named = true;
    result.chatId = chat.id;
  }
  workspace.activeId = chat.id;
  chat.snapshotOpen = true;
  chat.viewState = structuredClone(result.viewState);
  chat.activeStep = result.activeStep;
  if (!chat.messages.some(message => message.kind === 'snapshot')) chat.messages.push({kind:'snapshot', text:`Open ${result.name}.`});
  const folder = workspace.folders.find(item => item.id === chat.folderId);
  if (folder) folder.collapsed = false;
  restoreChat(); persistChats();
  $("#result").scrollIntoView({block:'start', behavior:'smooth'});
}
function startChat(folderId = null) {
  rememberChat();
  newChat(workspace, folderId);
  const folder = workspace.folders.find(f => f.id === folderId);
  if (folder) folder.collapsed = false;
  restoreChat();
  persistChats();
  $("#prompt").focus();
}
function switchChat(chatId) {
  if (chatId === workspace.activeId || !workspace.chats.some(c => c.id === chatId)) return;
  rememberChat();
  workspace.activeId = chatId;
  restoreChat();
  persistChats();
}
function messageHTML(message) {
  if (message.kind === 'snapshot') return `<div class="snapshot-message"><div class="user-message">${escape(message.text)}</div><button class="inline-action" data-open>Open saved results</button></div>`;
  return `<div class="exchange"><div class="user-message">${escape(message.text)}</div><div class="assistant-heading"><span class="assistant-avatar">r</span><b>Rebind</b></div><p class="assistant-message">New disease analysis isn’t connected yet. ${isDemo() ? "Example results use synthetic data." : `A saved analysis for <strong>${escape(manifest.disease.name)}</strong> is available.`}</p><button class="inline-action" data-open>${isDemo() ? "Open example results" : "Open saved analysis"}</button></div>`;
}
function renderConversation() {
  $("#conversation").innerHTML = currentChat().messages.map(messageHTML).join('');
}
function restoreChat() {
  const chat = currentChat();
  const defaults = initialState(manifest);
  state = {...defaults, ...chat.viewState};
  if (!['candidate','binding_site'].includes(state.view)) state.view = defaults.view;
  if (!['small_molecule','peptide','biologic'].includes(state.modality)) state.modality = defaults.modality;
  if (!['all','known','novel'].includes(state.novelty)) state.novelty = 'all';
  activeStep = chat.activeStep;
  opened = chat.snapshotOpen;
  $("#prompt").value = chat.draft;
  $("#thread-title").textContent = chat.name;
  renderConversation();
  const started = chat.messages.length > 0 || opened;
  $("#welcome").hidden = started;
  $("#main").classList.toggle('has-conversation', started);
  $("#result").hidden = !opened;
  renderChatList();
  if (opened) {setResultHeading(); renderResults();}
}
function submitPrompt(event) {
  event.preventDefault();
  const prompt = $("#prompt").value.trim();
  if (!prompt) return;
  const chat = currentChat();
  if (!chat.named) {chat.name = prompt.slice(0, 80); chat.named = true;}
  chat.messages.push({kind: 'prompt', text: prompt});
  chat.draft = '';
  chat.snapshotOpen = false;
  opened = false;
  restoreChat();
  persistChats();
  $("#conversation").lastElementChild?.scrollIntoView({block:'start', behavior:'smooth'});
}
function setResultHeading() {
  $("#origin-banner").textContent = isDemo() ? "Synthetic example" : manifest.data_origin === "computed" ? "Saved analysis" : "Reference example";
  $("#result-title").textContent = isDemo() ? "Example results" : manifest.disease.name;
  $("#result-description").textContent = isDemo() ? "Synthetic molecules and scores for testing." : "";
  $("#result-description").hidden = !isDemo();
}
function openSnapshot() {
  const chat = currentChat();
  if (!opened) {
    chat.messages.push({kind:'snapshot', text: isDemo() ? 'Show me the example results.' : `Open the saved ${manifest.disease.name} analysis.`});
    if (!chat.named) {chat.name = isDemo() ? 'Example results' : manifest.disease.name.slice(0,80); chat.named = true;}
    state = initialState(manifest);
    activeStep = isReference() ? 'signature' : 'candidates';
  }
  opened = true;
  begin();
  renderConversation();
  renderChatList();
  $("#thread-title").textContent = chat.name;
  $("#result").hidden = false;
  setResultHeading();
  renderResults();
  $("#result").scrollIntoView({block:'start', behavior:'smooth'});
}
function openOrganizer(kind, id = null) {
  organizerTarget = {kind, id};
  const item = kind === 'chat' ? workspace.chats.find(c => c.id === id) : workspace.folders.find(f => f.id === id);
  $("#organizer-title").textContent = kind === 'new-folder' ? 'New folder' : kind === 'chat' ? 'Edit chat' : 'Edit folder';
  $("#organizer-name").value = item?.name ?? '';
  $("#folder-field").hidden = kind !== 'chat';
  $("#organizer-folder").innerHTML = '<option value="">Unfiled</option>' + workspace.folders.map(f => `<option value="${f.id}">${escape(f.name)}</option>`).join('');
  $("#organizer-folder").value = item?.folderId ?? '';
  $("#organizer-delete").hidden = kind === 'new-folder';
  $("#organizer-delete").textContent = kind === 'chat' ? 'Delete chat' : 'Delete folder';
  $("#organizer-note").textContent = kind === 'folder' ? 'Deleting a folder keeps its chats in Unfiled.' : '';
  $("#organizer-error").textContent = '';
  $("#organizer").showModal();
  $("#organizer-name").focus();
}
function saveOrganizer(event) {
  event.preventDefault();
  const name = $("#organizer-name").value.trim();
  if (!name) {$("#organizer-error").textContent = 'Enter a name.'; return;}
  const {kind, id} = organizerTarget;
  if (kind === 'new-folder') addFolder(workspace, name);
  else if (kind === 'folder') workspace.folders.find(f => f.id === id).name = name.slice(0,80);
  else updateChat(workspace, id, name, $("#organizer-folder").value || null);
  persistChats(); renderChatList();
  $("#thread-title").textContent = currentChat().name;
  $("#organizer").close();
}
function deleteOrganizedItem() {
  const {kind, id} = organizerTarget;
  rememberChat();
  if (kind === 'folder') removeFolder(workspace, id);
  else removeChat(workspace, id);
  $("#organizer").close(); restoreChat(); persistChats();
}
function renderResults() {
  reconcile(manifest,state);
  state.chapter = 3;
  renderPipeline($("#pipeline-root"), manifest, activeStep, step => {activeStep = step; renderResults();});
  $("#result-panel").setAttribute("aria-labelledby", `step-${activeStep}`);
  $("#chapter-content").innerHTML = {candidates: board, literature: evidence, target: targetDetails, signature: bindingDetails}[activeStep]();
  $("#molecular-section").hidden = activeStep !== "signature";
  $("#validation-content").innerHTML = validation();
  $("#validation").hidden = !manifest.validation.length && !bundle.warnings.length;
  $("#modality")?.addEventListener("change", e=>{state.modality=e.target.value;renderResults();});
  $("#novelty")?.addEventListener("change", e=>{state.novelty=e.target.value;renderResults();});
  $("#rejected")?.addEventListener("change", e=>{state.showRejected=e.target.checked;renderResults();});
  if (activeStep === "signature") {
    if (!stage) stage = new MolecularStage($("#viewer-candidate"),$("#viewer-status"),key=>{state.residue=state.residue===key?null:key;renderResults();});
    const choices = rows(manifest, state);
    $("#structure-candidate").innerHTML = choices.length ? choices.map(c => `<option value="${escape(c.id)}" ${c.id === state.candidate ? "selected" : ""}>${escape(c.name)}</option>`).join("") : '<option value="">No candidate pose</option>';
    $("#structure-candidate").onchange = event => {state.candidate = event.target.value; renderResults();};
    $("#legend").innerHTML = '<span><i class="dot blue"></i>Candidate</span>' + (state.view === "binding_site" ? '<span><i class="dot teal"></i>Reference binding site</span><span>Target protein in grey</span>' : "");
    $("#residue-selection").textContent = state.residue ? `Selected residue: ${manifest.signature.residues.find(r=>r.key===state.residue)?.label ?? state.residue}` : "";
    document.querySelectorAll('[data-view]').forEach(el=>{el.classList.toggle('active',el.dataset.view===state.view);el.setAttribute('aria-pressed',String(el.dataset.view===state.view));});
    stage.draw(bundle,state).catch(()=>{$("#viewer-status").textContent="Structure unavailable. Saved scores are still accessible.";});
  }
  rememberChat();
  renderSavedResults();
}

function evidence() {
  const records = new Map();
  for (const e of manifest.evidence) {
    const key = e.url || e.id;
    if (!records.has(key)) records.set(key, {...e, claims: new Set()});
    records.get(key).claims.add(e.claim);
  }
  return `<div class="source-list">${[...records.values()].map(e => `<article><h3>${escape(e.title)}</h3>${[...e.claims].map(claim => `<p>${escape(claim)}</p>`).join("")}${e.url?.startsWith("https://") ? `<a href="${escape(e.url)}" target="_blank" rel="noopener noreferrer">${escape(e.id)} ↗</a>` : '<small>No citation supplied</small>'}</article>`).join("") || '<p>No sources supplied.</p>'}</div><details class="detail-disclosure"><summary>Run metadata</summary><dl class="provenance">${Object.entries(manifest.provenance).map(([key,value])=>`<dt>${escape(title(key))}</dt><dd>${escape(value)}</dd>`).join("")}</dl></details>`;
}

function targetDetails() {
  return `<div class="compact-summary"><h3>${escape(manifest.target.name)}</h3><small>${escape(manifest.target.id)} · ${escape(manifest.target.site_kind)}</small><p>${escape(manifest.target.rationale)}</p></div>`;
}

function bindingDetails() {
  const sig = manifest.signature;
  return `<details class="detail-disclosure binding-notes"><summary>${escape(shortSource[sig.source])} · ${sig.residues.length} residues</summary><p>${escape(sig.reason)}</p>${sig.source === "boltzgen_consensus" ? `<p>${sig.n_surviving} of ${sig.n_generated} designs retained; core threshold ${Math.round(sig.threshold * 100)}%.</p>` : ""}<p>Addressability: ${escape(title(sig.addressability.status))}. ${escape(sig.addressability.reason)}</p></details>`;
}

function candidateDetails(c) {
  if (!c) return "";
  const extras = `<p>${escape(c.approval)}</p>${!isReference() ? `<p>${escape(c.confidence.name)}: ${escape(c.confidence.value ?? "Not evaluated")} · ${escape(c.confidence.scale)}. ${escape(c.confidence.rule)}</p>` : ""}${c.caveats.length ? `<p>${c.caveats.map(escape).join(" ")}</p>` : ""}${c.decoys.status === "available" ? `<p>Decoy tie rule: ${escape(c.decoys.tie_rule)}. Percentile is not a binding probability.</p>` : `<p>${escape(c.decoys.reason)}</p>`}`;
  return `<details class="detail-disclosure"><summary>Selected candidate details</summary>${extras}</details>`;
}

function filters() {
  return `<div class="filters"><label>Modality <select id="modality">${["small_molecule", "peptide", "biologic"].map(value => `<option value="${value}" ${state.modality === value ? "selected" : ""}>${title(value)}</option>`).join("")}</select></label><label>Novelty <select id="novelty">${[["all", "All pairings"], ["novel", "Novel only"], ["known", "Known only"]].map(([value, text]) => `<option value="${value}" ${state.novelty === value ? "selected" : ""}>${text}</option>`).join("")}</select></label><label><input type="checkbox" id="rejected" ${state.showRejected ? "checked" : ""}> Include rejected / unassessed</label></div>`;
}

function board() {
  const visible = rows(manifest, state);
  const isSuppressed = suppressed(manifest, state.modality);
  return `<div class="section-head">${filters()}</div>${isSuppressed ? `<p class="subtle">${isReference() ? "Reference structure only. Candidate ranking has not been evaluated." : "Ranking unavailable; see the binding-site assessment."}</p>` : ""}<div class="panel table-scroll"><table><thead><tr><th>Rank</th><th>Molecule / pairing</th><th>Core coverage</th><th>Confidence gate</th><th>Decoy percentile</th></tr></thead><tbody>${visible.map(c => `<tr class="${c.id === state.candidate ? "selected" : ""}"><td>${!isSuppressed && c.rank !== null ? String(c.rank).padStart(2, "0") : "—"}</td><td><button class="table-button" data-candidate="${escape(c.id)}" aria-pressed="${c.id === state.candidate}">${escape(c.name)}<span class="subtext">${escape(title(c.novelty))}</span></button></td><td>${escape(percentage(c.coverage))}${c.coverage.value !== null ? `<span class="bar"><span style="width:${c.coverage.value * 100}%"></span></span>` : ""}</td><td><span class="tag ${c.confidence.gate === "passed" ? "" : "warning"}">${escape(title(c.confidence.gate))}</span></td><td>${c.decoys.status === "available" ? `${escape(c.decoys.percentile)}<span class="subtext">n = ${c.decoys.n}; within modality</span>` : "<span class='muted'>Uncalibrated</span>"}</td></tr>`).join("")}</tbody></table>${!visible.length ? "<div class='empty'>No candidates match these filters.</div>" : ""}</div><p class="subtle muted">Showing ${visible.length} of ${manifest.candidates.filter(c => c.modality === state.modality).length} records in this modality. Ordering supplied by the scientific pipeline: ${escape(visible[0]?.ranking_metric ?? "not available")}. Novel pairing means absent from the cited annotations, not proven novelty.</p><details class="detail-disclosure"><summary>View residue coverage</summary>${heatmap(visible)}</details>${candidateDetails(candidate())}`;
}

function heatmap(visible) {
  const core = manifest.signature.residues.filter(r => r.core);
  if (!core.length) return "<div class='notice'>No evaluable core signature. Coverage cannot be interpreted.</div>";
  return `<div class="section-head section"><h2>Residue coverage</h2><small>Select a residue to locate it in 3D.</small></div><div class="panel table-scroll"><table class="heatmap"><thead><tr><th>Candidate</th>${core.map(r => `<th><button data-residue="${escape(r.key)}">${escape(r.label)}<span class="subtext">${Math.round(r.frequency * 100)}% frequency</span></button></th>`).join("")}</tr></thead><tbody>${visible.map(c => `<tr><td><button class="table-button" data-candidate="${escape(c.id)}">${escape(c.name)}</button></td>${core.map(r => {const kind = c.engaged.includes(r.key) ? "engaged" : c.missed.includes(r.key) ? "missed" : "unavailable"; return `<td><button class="heat-cell ${kind} ${state.residue === r.key ? "active" : ""}" data-residue="${escape(r.key)}" aria-label="${escape(c.name)} / ${escape(r.label)}: ${kind}">${kind === "engaged" ? "+" : kind === "missed" ? "−" : "?"}</button></td>`;}).join("")}</tr>`).join("")}</tbody></table></div>`;
}

function validation() {
  return `<div class="cards">${manifest.validation.map(v => `<article class="card"><h3>${escape(v.label)}</h3><small>${escape(title(v.status))} · n = ${v.n}</small>${v.value !== null ? `<div class="stat">${escape(v.value)}</div>` : ""}<p>${escape(v.reason)}</p><small>${escape(v.source)}</small></article>`).join("")}</div>${bundle.warnings.map(w=>`<p>${escape(w)}</p>`).join("")}`;
}


try {
  const response = await fetch("bundle.json");
  if (!response.ok) throw new Error("Could not load the saved analysis.");
  bundle = await response.json();
  manifest = bundle.manifest;
  if (manifest?.schema_version !== "1.0.0" || !Array.isArray(manifest.candidates) || !bundle.display) throw new Error("This saved analysis needs to be rebuilt.");
  storageKey = `autorepurpose.chats.v1:${manifest.run_id}`;
  try {storage = window.localStorage;} catch {storage = null;}
  const saved = readWorkspace(storage, storageKey);
  workspace = saved.workspace;
  storageWritable = saved.writable;
  state = initialState(manifest);
  shell();
  restoreChat();
  persistChats();
} catch (error) {
  $("#app").innerHTML = `<main class="loading"><h1>Preview unavailable</h1><p>${escape(error.message)}</p><p>Build a visualization export and open it through the local server.</p></main>`;
}
