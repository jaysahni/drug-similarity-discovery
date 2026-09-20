export const chapters = [
  ["The question", "Evidence & rationale", "Start with a question worth asking.", "Trace the disease hypothesis to its target and supporting evidence."],
  ["The target", "Structure & site", "Find the place to intervene.", "Inspect the target structure and the addressability of its interface."],
  ["The signature", "Consensus hotspots", "A shared interface, made visible.", "Core residues describe recurring contacts—not experimentally proven energetic hotspots."],
  ["The match", "Covered & missed", "Same target. Different molecules.", "Compare target-side coverage, with structural confidence kept in view."],
  ["The shortlist", "Ranked hypotheses", "Make the next experiment informed.", "Inspect upstream rankings within each modality, not across incomparable molecules."],
  ["The reality check", "Validation & caveats", "Show the limits as clearly as the signal.", "Null comparisons, negative findings, and provenance belong beside every hypothesis."]
];

export function rows(manifest, state) {
  return manifest.candidates.filter(c => c.modality === state.modality)
    .filter(c => state.novelty === "all" || (state.novelty === "novel" ? c.novelty === "novel_pairing" : ["known_moa", "known_offtarget"].includes(c.novelty)))
    .filter(c => state.showRejected || c.confidence.gate === "passed")
    .sort((a, b) => (a.rank ?? Infinity) - (b.rank ?? Infinity));
}

export function percentage(measurement) {
  return measurement.status === "available" && measurement.value !== null ? `${Math.round(measurement.value * 100)}%` : "Not evaluated";
}

export function suppressed(manifest, modality) {
  return modality === "small_molecule" && manifest.signature.addressability.status !== "addressable";
}

export function initialState(manifest) {
  const preferred = manifest.candidates.find(c => c.rank !== null) ?? manifest.candidates[0];
  return {chapter: 0, modality: preferred?.modality ?? "small_molecule", novelty: "all", showRejected: false, candidate: preferred?.id ?? null, residue: null, view: "context"};
}

export function reconcile(manifest, state) {
  const visible = rows(manifest, state);
  if (!visible.some(c => c.id === state.candidate)) state.candidate = visible[0]?.id ?? null;
  return state;
}
