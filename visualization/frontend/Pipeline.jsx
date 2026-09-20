import React from 'react';
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

// Steps come from the manifest for a saved analysis, or are passed in for a recorded worked example.
export function Pipeline({manifest, steps: givenSteps, activeStep, onNavigate, panelId = 'result-panel'}) {
  const steps = givenSteps ?? pipelineSteps(manifest);
  function navigateWithKeys(event, index) {
    if (!['ArrowLeft', 'ArrowRight', 'ArrowUp', 'ArrowDown', 'Home', 'End'].includes(event.key)) return;
    event.preventDefault();
    const next = event.key === 'Home' ? 0 : event.key === 'End' ? steps.length - 1
      : (index + (['ArrowLeft', 'ArrowUp'].includes(event.key) ? -1 : 1) + steps.length) % steps.length;
    onNavigate(steps[next].id);
    event.currentTarget.closest('ol').querySelectorAll('[role="tab"]')[next].focus();
  }
  return (
    <nav className="pipeline" aria-label="Analysis pipeline">
      <ol className="pipeline-flow" role="tablist" aria-label="Analysis stages">
        {steps.map((step, index) => (
          <li className={`pipeline-step${step.available ? "" : " unavailable"}`} key={step.id} role="presentation">
            <button type="button" id={`step-${step.id}`} data-step={step.id}
              className={`pipeline-node${activeStep === step.id ? ' selected' : ''}`}
              role="tab" aria-selected={activeStep === step.id} aria-controls={panelId}
              tabIndex={activeStep === step.id ? 0 : -1}
              onKeyDown={event => navigateWithKeys(event, index)} onClick={() => onNavigate(step.id)}>
              <span className="pipeline-marker" aria-hidden="true">{index + 1}</span>
              <span className="pipeline-title">{step.title}</span>
            </button>
            {index < steps.length - 1 && <span className="pipeline-edge" aria-hidden="true"><i/></span>}
          </li>
        ))}
      </ol>
    </nav>
  );
}

const roots = new WeakMap();
export function renderPipeline(element, manifest, activeStep, onNavigate, options = {}) {
  let root = roots.get(element);
  if (!root) {root = createRoot(element); roots.set(element, root);}
  root.render(<Pipeline manifest={manifest} steps={options.steps} panelId={options.panelId} activeStep={activeStep} onNavigate={onNavigate}/>);
}
