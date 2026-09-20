// Worked examples recorded in SUCCESSES.md, rendered the way the portal shows a run.
// Pure functions with no DOM access, so the module is unit-testable under node:test.
const escape = value => String(value ?? "").replace(/[&<>"']/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));

export const exampleStepIds = ["disease", "protein", "search", "check"];
export const exampleStepTitles = {disease: "Disease → target", protein: "Target protein", search: "Candidate search", check: "Structural check"};
const roleLabels = {known: "known answer", hit: "repurposing hit", control: "control", reference: "reference", candidate: "candidate", boundary: "boundary", artifact: "artifact"};
const statusLabels = {recorded: "Recorded", partial: "Partly recorded", not_run: "Not run"};

export function validateExamples(data) {
  if (!data || data.version !== 1 || !Array.isArray(data.examples) || !data.examples.length) throw new Error("Worked examples are missing or unreadable.");
  const ids = new Set();
  for (const example of data.examples) {
    if (typeof example.id !== "string" || !/^[a-z0-9-]{1,40}$/.test(example.id) || ids.has(example.id)) throw new Error("Worked example ids must be unique slugs.");
    if (!example.disease || !example.target?.name || !Array.isArray(example.aliases) || !example.reply) throw new Error(`Worked example ${example.id} is incomplete.`);
    for (const id of exampleStepIds) {
      if (!statusLabels[example.steps?.[id]?.status]) throw new Error(`Worked example ${example.id} is missing the ${id} step.`);
    }
    ids.add(example.id);
  }
  return data.examples;
}

export const findExample = (examples, id) => examples.find(example => example.id === id) ?? null;

const normalise = text => String(text ?? "").toLowerCase().replace(/[^a-z0-9\s]/g, " ").replace(/\s+/g, " ").trim();

// Whole-word matching against the disease, target and aliases; the longest matching term wins.
export function matchExample(examples, prompt) {
  const text = ` ${normalise(prompt)} `;
  let best = null;
  for (const example of examples) {
    for (const alias of [example.disease, example.target.name, example.target.gene, ...example.aliases]) {
      const term = normalise(alias);
      if (term && text.includes(` ${term} `) && (!best || term.length > best.length)) best = {example, length: term.length};
    }
  }
  return best?.example ?? null;
}

export function exampleSteps(example) {
  return exampleStepIds.map(id => {
    const step = example.steps[id];
    return {id, title: exampleStepTitles[id], available: step.status !== "not_run", status: statusLabels[step.status], value: "", detail: step.intro ?? ""};
  });
}

const stats = list => list?.length ? `<div class="cards">${list.map(s => `<article class="card"><h3>${escape(s.label)}</h3><div class="stat">${escape(s.value)}</div><p>${escape(s.note)}</p></article>`).join("")}</div>` : "";
const notes = list => (list ?? []).map(note => `<p class="example-note">${escape(note)}</p>`).join("");
const notice = text => text ? `<div class="notice">${escape(text)}</div>` : "";
const intro = text => text ? `<p class="example-intro">${escape(text)}</p>` : "";
const tag = role => role ? `<span class="tag role-${escape(role)}">${escape(roleLabels[role] ?? role)}</span>` : "";

function diseaseStep(example, step) {
  const selected = `<div class="compact-summary"><h3>Selected target: ${escape(example.target.name)}</h3><small>gene ${escape(example.target.gene)}</small><p>${escape(step.selected)}</p></div>`;
  const rejected = step.rejected ? `<div class="section-head section"><h2>${escape(step.rejected.heading)}</h2><small>${escape(step.rejected.note)}</small></div><div class="panel table-scroll"><table><thead><tr><th>Protein</th><th>Why</th></tr></thead><tbody>${step.rejected.items.map(name => `<tr><td><b>${escape(name)}</b></td><td class="wrap">${escape(step.rejected.reason)}</td></tr>`).join("")}</tbody></table></div>` : "";
  const structure = step.structure ? `<div class="compact-summary section"><h3>Structure selection</h3><p>${escape(step.structure.examined)} Protein Data Bank entries examined; ${escape(step.structure.rejected)} rejected. ${escape(step.structure.reason)}</p></div>` : "";
  return `${intro(step.intro)}${stats(step.stats)}${selected}${notice(step.flag)}${rejected}${structure}${notes(step.notes)}${notice(step.unrecorded)}`;
}

function proteinStep(example, step) {
  const sequence = step.sequence ? `<p>${escape(step.sequence.length)}, ${escape(step.sequence.region)}. ${escape(step.sequence.note)}</p><small>Source: ${escape(step.source)}</small>` : "";
  return `${intro(step.intro)}<div class="compact-summary"><h3>${escape(step.name)}</h3><small>gene ${escape(step.gene)}</small>${sequence}${step.note ? `<p>${escape(step.note)}</p>` : ""}</div>${notes(step.notes)}${notice(step.unrecorded)}`;
}

function searchStep(example, step) {
  const rows = (step.rows ?? []).map(row => `<tr class="${row.role === "hit" ? "selected" : ""}"><td>${escape(row.rank)}</td><td><b>${escape(row.drug)}</b></td><td class="wrap">${escape(row.approved_for)}</td><td>${tag(row.role)}</td><td class="wrap">${escape(row.note)}</td></tr>`).join("");
  return `${intro(step.intro)}${stats(step.stats)}<div class="section-head section"><h2>What the blind search returned</h2><small>Query: ${escape(step.query)} · Library: ${escape(step.library)}</small></div><div class="panel table-scroll"><table><thead><tr><th>Rank</th><th>Drug</th><th>Approved for</th><th>Role</th><th>Note</th></tr></thead><tbody>${rows}</tbody></table></div>${notes(step.notes)}${notice(step.unrecorded)}`;
}

function checkStep(example, step) {
  if (step.status === "not_run") return `${intro(step.intro)}${notice(step.unrecorded)}`;
  const overlap = value => value === null || value === undefined ? "—" : `${escape(Number(value).toFixed(3))}<span class="bar"><span style="width:${Math.round(Number(value) * 100)}%"></span></span>`;
  const rows = (step.rows ?? []).map(row => `<tr class="${row.role === "hit" ? "selected" : ""}"><td><b>${escape(row.molecule)}</b>${row.detail ? `<span class="subtext">${escape(row.detail)}</span>` : ""}</td><td>${tag(row.role)}</td><td>${escape(row.confidence)}</td><td>${escape(row.contacts)}</td><td>${escape(row.shared)}</td><td>${overlap(row.overlap)}</td></tr>`).join("");
  return `${intro(step.intro)}<div class="section-head"><h2>Co-folded poses</h2><small>Reference: ${escape(step.reference)}</small></div><div class="panel table-scroll"><table><thead><tr><th>Molecule</th><th>Role</th><th>Confidence</th><th>Contacts</th><th>Shared with ${escape(step.reference)}</th><th>Overlap</th></tr></thead><tbody>${rows}</tbody></table></div>${notes(step.notes)}${notice(step.unrecorded)}`;
}

export function renderExampleStep(example, stepId) {
  const step = example.steps[stepId];
  if (!step) return "";
  return {disease: diseaseStep, protein: proteinStep, search: searchStep, check: checkStep}[stepId](example, step);
}

export function renderExampleDetails(example) {
  const list = (heading, items) => items?.length ? `<h3>${escape(heading)}</h3><ul>${items.map(item => `<li>${escape(item)}</li>`).join("")}</ul>` : "";
  const sources = example.sources?.length ? `<h3>Where the numbers live</h3><dl class="provenance">${example.sources.map(s => `<dt>${escape(s.label)}</dt><dd>${escape(s.path)}</dd>`).join("")}</dl>` : "";
  return `<p>Transcribed from SUCCESSES.md. Every number was computed by a script in this repository; the portal recomputes nothing.</p>${list("Why this matters", example.strengths)}${list("Limits", example.limits)}${sources}`;
}

export function renderExampleChips(examples, className = "example-chip") {
  return examples.map(example => `<button type="button" class="${escape(className)}" data-run-example="${escape(example.id)}">${escape(example.disease)}<span>${escape(example.target.name)}</span></button>`).join("");
}

// --- Co-folded structures for a recorded example -------------------------------
// The viewer only displays what the upstream contact extraction already decided.
export function validateStructures(data) {
  if (!data || data.version !== 1 || !data.examples) throw new Error("Structure data is unreadable.");
  for (const [id, entry] of Object.entries(data.examples)) {
    if (!Array.isArray(entry.molecules) || !entry.molecules.length) throw new Error(`No molecules for ${id}.`);
    for (const m of entry.molecules) {
      if (!m.id || typeof m.file !== "string" || !/^poses\/[A-Za-z0-9_.-]+\.pdb$/.test(m.file)) throw new Error(`Bad pose path for ${m.id}.`);
      if (!Array.isArray(m.ligand_atoms) || !m.ligand_atoms.length) throw new Error(`No ligand atoms for ${m.id}.`);
    }
  }
  return data;
}

export const structuresFor = (data, exampleId) => data?.examples?.[exampleId] ?? null;
// With nothing chosen yet, open on the hit: it is the claim the panel exists to show.
export const moleculeIn = (entry, id) => entry?.molecules.find(m => m.id === id)
  ?? entry?.molecules.find(m => m.role === "hit") ?? entry?.molecules[0] ?? null;

// Shape the pose as the MolecularStage already understands, so the tested viewer is reused unchanged.
export function structureBundle(molecule, pdb, colors) {
  return {
    manifest: {
      candidates: [{id: molecule.id, structure_id: molecule.id}],
      signature: {residues: Object.keys(molecule.contact_atoms ?? {}).map(key => ({key, label: key}))},
    },
    display: {[molecule.id]: {pdb, ligand_atoms: molecule.ligand_atoms, residue_atoms: molecule.contact_atoms ?? {}}},
    colors,
  };
}

export function renderStructurePanel(entry, selectedId, view) {
  if (!entry) return "";
  const molecule = moleculeIn(entry, selectedId);
  const option = m => `<option value="${escape(m.id)}" ${m.id === molecule.id ? "selected" : ""}>${escape(m.name)} · ${escape(roleLabels[m.role] ?? m.role)}</option>`;
  const views = [["candidate", "Molecule only"], ["binding_site", "Molecule in the pocket"]];
  const stat = (label, value) => `<div class="pose-stat"><span>${escape(label)}</span><b>${escape(value)}</b></div>`;
  const figures = [
    stat("Confidence", molecule.confidence),
    stat("Contact residues", molecule.contacts),
    molecule.shared === null ? stat("Role", "reference pose") : stat(`Shared with ${entry.reference}`, `${molecule.shared} of ${molecule.contacts}`),
    molecule.overlap === null ? "" : stat("Overlap", molecule.overlap.toFixed(3)),
  ].join("");
  return `<div class="panel pose-panel"><div class="panel-heading"><h2>Co-folded structure</h2><div class="segmented" role="group" aria-label="Structure view">${views.map(([id, label]) => `<button data-example-view="${id}" class="${view === id ? "active" : ""}" aria-pressed="${view === id}">${label}</button>`).join("")}</div></div>
    <label class="molecular-label" for="example-molecule">Molecule <select id="example-molecule">${entry.molecules.map(option).join("")}</select></label>
    <div id="example-viewer" class="viewer" role="img" aria-label="Co-folded structure of ${escape(molecule.name)} with ${escape(entry.target)}"></div>
    <div class="pose-stats">${figures}</div>
    <div class="legend"><span><i class="dot blue"></i>${escape(molecule.name)}</span>${view === "binding_site" ? `<span><i class="dot teal"></i>Residues it touches</span><span>${escape(entry.target)} in grey</span>` : ""}</div>
    ${entry.alignment?.status === "failed" ? `<p class="pose-caveat">${escape(entry.alignment.reason)}</p>` : ""}
    <div id="example-viewer-status" class="viewer-status" role="status"></div></div>
    <p class="example-note">${escape(entry.caption)}</p>`;
}
