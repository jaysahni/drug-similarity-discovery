import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {test} from 'node:test';

const source = await readFile(new URL('../web/viewer.js', import.meta.url), 'utf8');
const {MolecularStage} = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
const bundle = {
  manifest: {candidates: [{id: 'drug', structure_id: 'pose'}], signature: {residues: [{key: 'A:1', label: 'A:1'}]}},
  display: {pose: {pdb: 'pose', ligand_atoms: [10, 11], residue_atoms: {'A:1': [1, 2]}}},
  colors: {ligand: 'blue', target: 'grey', engaged: 'teal', selected: 'purple'},
};
function stage() {
  const calls = [];
  const value = Object.create(MolecularStage.prototype);
  value.ready = true;
  value.revision = 0;
  value.status = {};
  value.onResidue = key => calls.push(['residue', key]);
  value.viewer = Object.fromEntries(['removeAllSurfaces', 'removeAllLabels', 'removeAllModels', 'addModel', 'setStyle', 'addStyle', 'addLabel', 'setClickable', 'setProjection', 'zoomTo', 'render', 'addSurface'].map(name => [name, (...args) => {calls.push([name, ...args]);}]));
  value.viewer.selectedAtoms = () => [{serial: 1}];
  return {value, calls};
}
globalThis.window = {$3Dmol: {SurfaceType: {VDW: 'VDW'}}};

test('candidate mode hides the target and frames only candidate atoms', async () => {
  const {value, calls} = stage();
  await value.draw(bundle, {candidate: 'drug', view: 'candidate', residue: 'A:1'});
  assert.deepEqual(calls.filter(c => c[0] === 'setStyle'), [
    ['setStyle', {}, {}],
    ['setStyle', {serial: [10, 11]}, {stick: {color: 'blue', radius: 0.24}}],
  ]);
  assert.deepEqual(calls.find(c => c[0] === 'zoomTo')[1], {serial: [10, 11]});
  assert.equal(calls.some(c => c[0] === 'addSurface' || c[0] === 'addLabel'), false);
});

test('binding-site mode uses the same pose, maps residues, and frames both selections', async () => {
  const {value, calls} = stage();
  await value.draw(bundle, {candidate: 'drug', view: 'binding_site', residue: 'A:1'});
  assert.deepEqual(calls.filter(c => c[0] === 'addModel'), [['addModel', 'pose', 'pdb']]);
  assert.deepEqual(calls.find(c => c[0] === 'zoomTo')[1], {serial: [10, 11, 1, 2]});
  assert.deepEqual(calls.find(c => c[0] === 'addSurface')[3], {serial: [1, 2]});
  assert.equal(calls.find(c => c[0] === 'addLabel')[1], 'A:1');
  calls.find(c => c[0] === 'setClickable')[3]({serial: 1});
  assert.deepEqual(calls.at(-1), ['residue', 'A:1']);
  await value.draw(bundle, {candidate: 'drug', view: 'candidate'});
  assert.equal(calls.filter(c => c[0] === 'removeAllModels').length, 2);
  assert.equal(value.status.textContent, 'Candidate molecule only.');
});

test('missing candidate geometry does not substitute a reference protein', async () => {
  for (const display of [{}, {pose: {...bundle.display.pose, ligand_atoms: []}}]) {
    const {value, calls} = stage();
    await value.draw({...bundle, display}, {candidate: 'drug', view: 'binding_site'});
    assert.equal(calls.some(c => c[0] === 'addModel'), false);
    assert.ok(value.status.textContent);
  }
});
