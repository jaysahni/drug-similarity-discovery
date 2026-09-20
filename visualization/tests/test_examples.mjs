import assert from 'node:assert/strict';
import {existsSync} from 'node:fs';
import {readFile} from 'node:fs/promises';
import {pathToFileURL} from 'node:url';
import {test} from 'node:test';

const source = await readFile(new URL('../web/examples.js', import.meta.url), 'utf8');
const {validateExamples, matchExample, exampleSteps, exampleStepIds, renderExampleStep, renderExampleDetails, renderExampleChips} = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
const data = JSON.parse(await readFile(new URL('../web/examples.json', import.meta.url), 'utf8'));
const examples = validateExamples(data);

test('the shipped worked examples validate and cover every pipeline step', () => {
  assert.deepEqual(examples.map(e => e.id), ['kdr', 'egfr', 'padi4']);
  for (const example of examples) for (const id of exampleStepIds) assert.ok(example.steps[id].status, `${example.id} ${id}`);
  assert.throws(() => validateExamples({version: 1, examples: [{id: 'x'}]}));
  assert.throws(() => validateExamples({version: 1, examples: [{...examples[0], id: 'bad id'}]}));
  assert.throws(() => validateExamples({version: 1, examples: [examples[0], examples[0]]}));
  assert.throws(() => validateExamples({version: 2, examples}));
});

test('prompts match examples by whole words, longest term first', () => {
  assert.equal(matchExample(examples, 'Colorectal cancer')?.id, 'kdr');
  assert.equal(matchExample(examples, 'which approved drugs might treat lung cancer?')?.id, 'egfr');
  assert.equal(matchExample(examples, 'RA')?.id, 'padi4');
  assert.equal(matchExample(examples, 'VEGFR-2 inhibitors')?.id, 'kdr');
  assert.equal(matchExample(examples, 'brainstorm'), null);
  assert.equal(matchExample(examples, 'Diabetes'), null);
  assert.equal(matchExample(examples, 'Open example results'), null);
  assert.equal(matchExample(examples, ''), null);
});

test('steps that were not run are marked unavailable in the pipeline', () => {
  assert.deepEqual(exampleSteps(examples[0]).map(s => s.available), [true, true, true, true]);
  assert.deepEqual(exampleSteps(examples[1]).map(s => s.available), [true, true, true, false]);
  assert.deepEqual(exampleSteps(examples[2]).map(s => s.available), [true, true, true, false]);
  assert.deepEqual(exampleSteps(examples[0]).map(s => s.id), exampleStepIds);
});

test('rendering carries the recorded numbers and escapes markup', () => {
  const kdr = examples[0];
  assert.match(renderExampleStep(kdr, 'disease'), /16,299/);
  assert.match(renderExampleStep(kdr, 'disease'), /MSH2/);
  assert.match(renderExampleStep(kdr, 'protein'), /329 residues/);
  assert.match(renderExampleStep(kdr, 'search'), /mebendazole/);
  assert.match(renderExampleStep(kdr, 'check'), /0\.950/);
  assert.match(renderExampleStep(examples[1], 'check'), /not run/i);
  assert.match(renderExampleStep(examples[2], 'search'), /pentamidine/);
  assert.match(renderExampleDetails(kdr), /SUCCESSES\.md/);
  assert.equal(renderExampleStep(kdr, 'nowhere'), '');
  const hostile = structuredClone(kdr);
  hostile.disease = '<img src=x onerror=alert(1)>';
  hostile.steps.search.rows[0].drug = '<b>x</b>';
  hostile.steps.search.stats[0].value = '<i>1</i>';
  const html = renderExampleStep(hostile, 'search') + renderExampleChips([hostile]);
  assert.equal(html.includes('<img'), false);
  assert.equal(html.includes('<b>x'), false);
  assert.equal(html.includes('<i>1'), false);
});

const successes = process.env.SUCCESSES_MD ? pathToFileURL(process.env.SUCCESSES_MD) : new URL('../../SUCCESSES.md', import.meta.url);
test('every recorded number appears verbatim in SUCCESSES.md', {skip: existsSync(successes) ? false : 'SUCCESSES.md is not checked out'}, async () => {
  const doc = await readFile(successes, 'utf8');
  const numeric = value => typeof value === 'string' && /^\d[\d,.]*×?$/.test(value);
  const missing = [];
  for (const example of examples) {
    for (const step of Object.values(example.steps)) {
      for (const stat of step.stats ?? []) if (numeric(stat.value) && !doc.includes(stat.value)) missing.push(`${example.id}: ${stat.value}`);
      for (const row of step.rows ?? []) {
        for (const field of ['rank', 'confidence', 'contacts', 'shared']) if (numeric(row[field]) && !doc.includes(row[field])) missing.push(`${example.id}: ${row[field]}`);
        if (typeof row.overlap === 'number' && !doc.includes(row.overlap.toFixed(3))) missing.push(`${example.id}: overlap ${row.overlap}`);
        if (row.drug && !row.drug.startsWith('Known') && !doc.toLowerCase().includes(row.drug.split(',')[0].toLowerCase())) missing.push(`${example.id}: ${row.drug}`);
        if (row.molecule && !doc.toLowerCase().includes(row.molecule.toLowerCase())) missing.push(`${example.id}: ${row.molecule}`);
      }
    }
  }
  assert.deepEqual(missing, []);
});

const {validateStructures, structuresFor, moleculeIn, structureBundle, renderStructurePanel} = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
const structures = validateStructures(JSON.parse(await readFile(new URL('../web/example-structures.json', import.meta.url), 'utf8')));

test('the shipped co-folded poses validate and match the documented table', () => {
  const entry = structuresFor(structures, 'kdr');
  assert.equal(entry.reference, 'axitinib');
  assert.deepEqual(entry.molecules.map(m => m.id), ['axitinib', 'mebendazole', 'niclosamide', 'paracetamol', 'warfarin', 'furosemide']);
  const expected = {
    axitinib:    {confidence: 0.992, contacts: 20, shared: null, overlap: null},
    mebendazole: {confidence: 0.99,  contacts: 19, shared: 19,   overlap: 0.95},
    niclosamide: {confidence: 0.978, contacts: 17, shared: 11,   overlap: 0.423},
    paracetamol: {confidence: 0.97,  contacts: 11, shared: 9,    overlap: 0.409},
    warfarin:    {confidence: 0.953, contacts: 17, shared: 11,   overlap: 0.423},
    furosemide:  {confidence: 0.883, contacts: 16, shared: 14,   overlap: 0.636},
  };
  for (const m of entry.molecules) {
    assert.deepEqual({confidence: m.confidence, contacts: m.contacts, shared: m.shared, overlap: m.overlap}, expected[m.id], m.id);
    // Every contact residue carries atoms to highlight, and the ligand is never empty.
    assert.equal(Object.keys(m.contact_atoms).length, m.contacts, `${m.id} contact atoms`);
    assert.ok(m.ligand_atoms.length > 0);
    assert.ok(Object.values(m.contact_atoms).every(list => Array.isArray(list) && list.length));
  }
});

test('structure data is rejected when a pose path or ligand is unusable', () => {
  const entry = structuresFor(structures, 'kdr');
  const broken = path => ({version: 1, examples: {kdr: {...entry, molecules: [{...entry.molecules[0], ...path}]}}});
  assert.throws(() => validateStructures(broken({file: '../../etc/passwd'})));
  assert.throws(() => validateStructures(broken({file: 'poses/x.txt'})));
  assert.throws(() => validateStructures(broken({ligand_atoms: []})));
  assert.throws(() => validateStructures({version: 2, examples: {}}));
  assert.equal(structuresFor(structures, 'egfr'), null);
});

test('the pose is shaped for the existing viewer, and the panel reflects the selection', () => {
  const entry = structuresFor(structures, 'kdr');
  const meb = moleculeIn(entry, 'mebendazole');
  const shaped = structureBundle(meb, 'PDBTEXT', {ligand: 'blue'});
  assert.equal(shaped.display.mebendazole.pdb, 'PDBTEXT');
  assert.deepEqual(shaped.manifest.candidates, [{id: 'mebendazole', structure_id: 'mebendazole'}]);
  assert.equal(shaped.manifest.signature.residues.length, meb.contacts);
  assert.equal(shaped.display.mebendazole.ligand_atoms, meb.ligand_atoms);
  assert.equal(moleculeIn(entry, 'nonexistent').id, 'mebendazole', 'opens on the hit when nothing is chosen');
  assert.equal(moleculeIn(entry, null).role, 'hit');
  const html = renderStructurePanel(entry, 'warfarin', 'binding_site');
  assert.match(html, /warfarin/);
  assert.match(html, /0\.423/);
  assert.match(html, /aria-pressed="true"[^>]*>Molecule in the pocket|Molecule in the pocket/);
  assert.equal(renderStructurePanel(null, 'x', 'candidate'), '');
  // The poses failed the alignment gate, so the panel must say so rather than invite a cross-molecule comparison.
  assert.equal(entry.alignment.status, 'failed');
  assert.match(html, /own frame|synchronized comparison stays disabled/);
  assert.equal(/land somewhere else/.test(JSON.stringify(entry)), false);
});

test('the structure panel is not specific to one example', () => {
  // A second example id must work with no code change; only data decides what can be shown.
  const kdr = structuresFor(structures, 'kdr');
  const second = {
    version: 1,
    examples: {
      kdr,
      other: {
        reference: 'ADP', target: 'MSH2', step: 'check', caption: 'A different target entirely.',
        molecules: [
          {id: 'adp', name: 'ADP', role: 'reference', detail: '', file: 'poses/OTHER_adp.pdb',
           confidence: 0.972, contacts: 16, shared: null, overlap: null,
           ligand_atoms: [1, 2], contact_atoms: {'A:5': [3]}},
          {id: 'warfarin', name: 'warfarin', role: 'control', detail: 'blood thinner', file: 'poses/OTHER_warfarin.pdb',
           confidence: 0.547, contacts: 17, shared: 2, overlap: 0.065,
           ligand_atoms: [4], contact_atoms: {'A:9': [5]}},
        ],
      },
    },
  };
  assert.doesNotThrow(() => validateStructures(second));
  const entry = structuresFor(validateStructures(second), 'other');
  assert.equal(entry.target, 'MSH2');
  // With no hit in the set it opens on the first molecule, and the panel renders that target's numbers.
  assert.equal(moleculeIn(entry, null).id, 'adp');
  const html = renderStructurePanel(entry, 'warfarin', 'binding_site');
  assert.match(html, /MSH2/);
  assert.match(html, /0\.065/);
  assert.match(html, /Shared with ADP/);
  assert.equal(structureBundle(moleculeIn(entry, 'warfarin'), 'TEXT', {}).display.warfarin.pdb, 'TEXT');
  // Still keyed independently: the colorectal entry is untouched by the new one.
  assert.equal(structuresFor(validateStructures(second), 'kdr').molecules.length, 6);
});
