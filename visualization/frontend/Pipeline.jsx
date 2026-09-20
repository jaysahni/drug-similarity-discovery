import React, {useState} from 'react';
import {createRoot} from 'react-dom/client';

const sourceNames = {
  boltzgen_consensus: 'BoltzGen consensus',
  p2rank_geometry: 'Pocket geometry',
  known_ligand: 'Known-ligand reference',
  hybrid: 'Hybrid signature',
};
const countLabel = (n, singular, plural = `${singular}s`) => `${n} ${n === 1 ? singular : plural}`;

export function pipelineSteps(manifest) {
  const sources = manifest.evidence.filter(e => e.url?.startsWith('https://')).length;
  const hasTarget = Boolean(manifest.target.structure_id);
  const hasSite = manifest.signature.status === 'available';
  const reference = manifest.data_origin === 'reference_example';
  const ranked = manifest.candidates.filter(c => c.rank !== null).length;
  return [
    {
      id: 'literature', title: 'Literature', icon: 'papers',
      value: countLabel(sources, 'linked source'),
      status: sources ? 'Available' : 'Not supplied', available: sources > 0,
      detail: sources ? 'Sources supplied with this analysis. No live literature search was run.' : 'This analysis has no linked sources.',
      action: 'View sources', tab: 'sources', connection: 'Evidence',
    },
    {
      id: 'target', title: 'Target protein', icon: 'target',
      value: manifest.target.name,
      status: hasTarget ? 'Structure supplied' : 'No structure', available: hasTarget,
      detail: manifest.target.rationale,
      action: 'View structure', tab: 'structure', view: 'binding_site', connection: 'Binding site',
    },
    {
      id: 'signature', title: 'Binding site', icon: 'site',
      value: sourceNames[manifest.signature.source] ?? 'Not supplied',
      status: hasSite ? countLabel(manifest.signature.residues.length, 'residue') : 'Not evaluated', available: hasSite,
      detail: manifest.signature.reason,
      action: 'View binding site', tab: 'structure', view: 'binding_site', connection: 'Comparison',
    },
    {
      id: 'candidates', title: 'Candidates', icon: 'molecule',
      value: countLabel(manifest.candidates.length, 'record'),
      status: reference ? 'Reference only' : ranked ? `${ranked} ranked` : 'Unranked', available: manifest.candidates.length > 0,
      detail: reference ? 'The recorded ligand is shown as a structural reference. No candidate ranking was performed.' : `${countLabel(ranked, 'candidate')} ranked in the saved results.`,
      action: 'View candidates', tab: 'candidates',
    },
  ];
}

export function Pipeline({manifest, onNavigate}) {
  const [selected, setSelected] = useState(null);
  const steps = pipelineSteps(manifest);
  const current = steps.find(step => step.id === selected);
  const origin = manifest.data_origin === 'reference_example' ? 'Reference example' : manifest.data_origin === 'synthetic_fixture' ? 'Synthetic example' : 'Saved analysis';
  return (
    <section className="pipeline" aria-label="Analysis pipeline">
      <header className="pipeline-heading"><h3>Pipeline</h3><span>{origin}</span></header>
      <ol className="pipeline-flow">
        {steps.map((step, index) => (
          <li className={`pipeline-step${step.available ? "" : " unavailable"}`} key={step.id}>
            <button type="button" className={`pipeline-node${selected === step.id ? ' selected' : ''}`}
              aria-expanded={selected === step.id} aria-controls="pipeline-detail"
              onClick={() => setSelected(selected === step.id ? null : step.id)}>
              <span className="pipeline-marker" aria-hidden="true">{index + 1}</span>
              <span className="pipeline-title">{step.title}</span>
              <span className="pipeline-value">{step.value}</span>
            </button>
            {index < steps.length - 1 && <span className="pipeline-edge" aria-hidden="true"><i/></span>}
          </li>
        ))}
      </ol>
      <div id="pipeline-detail" className="pipeline-detail" hidden={!current}>
        {current && <><div><div className="pipeline-detail-heading"><strong>{current.title}</strong><span>{current.status}</span></div><p>{current.detail}</p></div><button type="button" onClick={() => onNavigate(current.tab, current.view)}>{current.action} <span aria-hidden="true">→</span></button></>}
      </div>
    </section>
  );
}

const roots = new WeakMap();
export function renderPipeline(element, manifest, onNavigate) {
  let root = roots.get(element);
  if (!root) {root = createRoot(element); roots.set(element, root);}
  root.render(<Pipeline manifest={manifest} onNavigate={onNavigate}/>);
}
