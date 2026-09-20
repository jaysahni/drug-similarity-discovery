export class MolecularStage {
  constructor(left, right, status, onResidue) {
    this.elements = [left, right];
    this.status = status;
    this.onResidue = onResidue;
    this.viewers = [];
    this.revision = 0;
    this.ready = false;
    try {
      if (!window.$3Dmol) throw new Error("3Dmol asset unavailable");
      const probe = document.createElement("canvas");
      const gl = probe.getContext("webgl2") || probe.getContext("webgl");
      if (!gl) throw new Error("WebGL unavailable");
      gl.getExtension("WEBGL_lose_context")?.loseContext();
      this.viewers = this.elements.map(element => window.$3Dmol.createViewer(element, {backgroundColor: "white", antialias: true}));
      this.viewers.forEach((v, i) => this.elements[i].querySelector("canvas")?.addEventListener("webglcontextlost", event => {
        event.preventDefault();
        this.fail("The graphics context was lost. Tables and residue explanations remain available; reload to restore 3D.");
      }));
      this.elements.forEach((element, index) => {
        for (const eventName of ["mousemove", "mouseup", "wheel", "touchmove"]) element.addEventListener(eventName, () => {
          if (this.synchronized) requestAnimationFrame(() => {
            if (!this.synchronized || !this.ready) return;
            this.viewers[1 - index].setView(this.viewers[index].getView());
            this.viewers[1 - index].render();
          });
        }, {passive: true});
      });
      this.resizeObserver = new ResizeObserver(() => this.viewers.forEach(v => {v.resize(); v.render();}));
      this.elements.forEach(el => this.resizeObserver.observe(el));
      this.ready = true;
    } catch (error) {
      this.fail("3D unavailable. Molecular rendering needs WebGL and the bundled 3Dmol asset. All evidence and results below remain usable.");
    }
  }

  fail(message) {
    this.ready = false;
    this.status.textContent = message;
    this.elements.forEach(el => {
      el.replaceChildren();
      const text = document.createElement("p");
      text.className = "viewer-empty";
      text.textContent = "Structure view unavailable";
      el.append(text);
    });
  }

  async draw(bundle, state) {
    if (!this.ready) return;
    const revision = ++this.revision;
    const manifest = bundle.manifest;
    const candidate = manifest.candidates.find(c => c.id === state.candidate);
    const reference = state.chapter >= 2 && manifest.design_reference ? manifest.design_reference.structure_id : manifest.target.structure_id;
    const assets = [bundle.display[reference], bundle.display[candidate?.structure_id]];
    const sameFrame = assets[0] && assets[1] && assets[0].frame_id === assets[1].frame_id;
    this.status.textContent = sameFrame ? "Target and candidate share a declared coordinate frame. Views are linked." : "Independent coordinate frames or missing pose. No overlay or synchronized alignment is implied.";
    const jobs = this.viewers.map(async (viewer, index) => {
      viewer.removeAllSurfaces();
      viewer.removeAllLabels();
      viewer.removeAllModels();
      const asset = assets[index];
      if (!asset) {
        viewer.render();
        return;
      }
      viewer.addModel(asset.pdb, "pdb");
      viewer.setStyle({}, {cartoon: {color: bundle.colors.target}, line: {color: bundle.colors.target, opacity: 0.4}});
      const serialToKey = new Map();
      const site = [];
      for (const residue of manifest.signature.residues) {
        const serial = asset.residue_atoms[residue.key] ?? [];
        serial.forEach(n => serialToKey.set(n, residue.key));
        site.push(...serial);
        let color = bundle.colors.target;
        if (state.view === "frequency" || index === 0 && state.chapter >= 2) {
          const f = residue.frequency;
          color = `rgb(${Math.round(235 * (1 - f))},${Math.round(245 - 118 * f)},${Math.round(245 - 130 * f)})`;
        } else if (index === 1 && candidate && ["coverage", "missed"].includes(state.view)) {
          color = candidate.engaged.includes(residue.key) ? bundle.colors.engaged : candidate.missed.includes(residue.key) ? bundle.colors.missed : bundle.colors.target;
        } else if (state.view === "site") color = bundle.colors.engaged;
        if (serial.length && state.view !== "context") {
          viewer.setStyle({serial}, {cartoon: {color}, stick: {color, radius: residue.key === state.residue ? 0.3 : 0.15}});
        }
        if (residue.key === state.residue && serial.length) {
          viewer.addStyle({serial}, {sphere: {color: bundle.colors.selected, radius: 0.32}});
          const atom = viewer.selectedAtoms({serial: [serial[0]]})[0];
          if (atom) viewer.addLabel(residue.label, {position: atom, backgroundColor: "white", fontColor: "#172b3a", borderColor: bundle.colors.selected, borderThickness: 1, fontSize: 13});
        }
      }
      if (asset.ligand_atoms.length) {
        const showDesign = index === 0 && state.chapter >= 2 && manifest.design_reference;
        viewer.setStyle({serial: asset.ligand_atoms}, index === 1 ? {stick: {color: bundle.colors.ligand, radius: 0.24}} : showDesign ? {cartoon: {color: bundle.colors.design}, stick: {color: bundle.colors.design, radius: 0.1}} : {});
      }
      viewer.setClickable({}, true, atom => {
        const key = serialToKey.get(atom.serial);
        if (key) this.onResidue(key);
      });
      viewer.setProjection("orthographic");
      viewer.zoomTo(state.view !== "context" && site.length ? {serial: site} : {});
      viewer.render();
      if (site.length && state.view !== "context") {
        try {
          await viewer.addSurface(window.$3Dmol.SurfaceType.VDW, {opacity: 0.15, color: bundle.colors.target}, {serial: site});
          if (revision === this.revision) viewer.render();
        } catch (error) {
          if (revision === this.revision) this.status.textContent = "Surface rendering unavailable; atom and residue views remain active.";
        }
      }
    });
    await Promise.all(jobs);
    if (revision !== this.revision) return;
    this.synchronized = Boolean(sameFrame);
    if (sameFrame) {
      this.viewers[1].setView(this.viewers[0].getView());
      this.viewers[1].render();
    }
  }
}
