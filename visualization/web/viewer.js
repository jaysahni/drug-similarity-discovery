export class MolecularStage {
  constructor(element, status, onResidue) {
    this.element = element;
    this.status = status;
    this.onResidue = onResidue;
    this.revision = 0;
    this.ready = false;
    try {
      if (!window.$3Dmol) throw new Error("3Dmol asset unavailable");
      const probe = document.createElement("canvas");
      const gl = probe.getContext("webgl2") || probe.getContext("webgl");
      if (!gl) throw new Error("WebGL unavailable");
      gl.getExtension("WEBGL_lose_context")?.loseContext();
      this.viewer = window.$3Dmol.createViewer(element, {backgroundColor: "white", antialias: true});
      element.querySelector("canvas")?.addEventListener("webglcontextlost", event => {
        event.preventDefault();
        this.fail("3D rendering stopped. Reload to restore it.");
      });
      this.resizeObserver = new ResizeObserver(() => {
        if (this.ready) {this.viewer.resize(); this.viewer.render();}
      });
      this.resizeObserver.observe(element);
      this.ready = true;
    } catch (error) {
      this.fail("3D unavailable. Saved scores are still accessible.");
    }
  }

  fail(message) {
    this.ready = false;
    this.status.textContent = message;
    this.element.replaceChildren();
    const text = document.createElement("p");
    text.className = "viewer-empty";
    text.textContent = "Structure view unavailable";
    this.element.append(text);
  }

  async draw(bundle, state) {
    if (!this.ready) return;
    const revision = ++this.revision;
    const viewer = this.viewer;
    const candidate = bundle.manifest.candidates.find(c => c.id === state.candidate);
    const asset = bundle.display[candidate?.structure_id];
    const withSite = state.view === "binding_site";
    viewer.removeAllSurfaces();
    viewer.removeAllLabels();
    viewer.removeAllModels();
    if (!asset || !asset.ligand_atoms.length) {
      this.status.textContent = !asset ? "No structure available for this candidate." : "Candidate atoms are not identified in this structure.";
      viewer.render();
      return;
    }
    viewer.addModel(asset.pdb, "pdb");
    // Hide every atom first; the candidate view contains only the supplied binder selection.
    viewer.setStyle({}, {});
    const site = [];
    const serialToKey = new Map();
    if (withSite) {
      viewer.setStyle({serial: asset.ligand_atoms, invert: true}, {cartoon: {color: bundle.colors.target}});
      for (const residue of bundle.manifest.signature.residues) {
        const serial = asset.residue_atoms[residue.key] ?? [];
        if (!serial.length) continue;
        site.push(...serial);
        serial.forEach(n => serialToKey.set(n, residue.key));
        viewer.setStyle({serial}, {cartoon: {color: bundle.colors.engaged}, stick: {color: bundle.colors.engaged, radius: 0.15}});
        if (residue.key === state.residue) {
          viewer.addStyle({serial}, {sphere: {color: bundle.colors.selected, radius: 0.32}});
          const atom = viewer.selectedAtoms({serial: [serial[0]]})[0];
          if (atom) viewer.addLabel(residue.label, {position: atom, backgroundColor: "white", fontColor: "#172b3a", fontSize: 13});
        }
      }
      viewer.setClickable({serial: site}, true, atom => {
        const key = serialToKey.get(atom.serial);
        if (key) this.onResidue(key);
      });
    }
    viewer.setStyle({serial: asset.ligand_atoms}, {stick: {color: bundle.colors.ligand, radius: 0.24}});
    viewer.setProjection("orthographic");
    viewer.zoomTo({serial: [...asset.ligand_atoms, ...site]});
    this.status.textContent = withSite ? (site.length ? "Binding site from the saved reference signature." : "No binding site mapped in this candidate structure.") : "Candidate molecule only.";
    viewer.render();
    if (site.length) {
      try {
        await viewer.addSurface(window.$3Dmol.SurfaceType.VDW, {opacity: 0.15, color: bundle.colors.target}, {serial: site});
        if (revision === this.revision) viewer.render();
      } catch (error) {
        if (revision === this.revision) this.status.textContent = "Surface unavailable; binding-site atoms are shown.";
      }
    }
  }
}
