"""Package the existing 3VHE complex as an unranked structural reference.

Run from the repository root:
    visualization/.venv/bin/python -m visualization.reference --out visualization/out/reference-source

Uses existing upstream contact annotations for display only; does not compute
contacts, scores, predictions, or modify scientific outputs.
"""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

from Bio.PDB import PDBParser
from .bundle import BundleError, validate

ROOT = Path(__file__).resolve().parent.parent


def make_reference(destination):
    source = ROOT / 'results/pipeline/colorectal-cancer/target'
    pdb = source / 'structure/3VHE_A_holo.pdb'
    contacts = json.loads((source / 'known_ligand_contacts.json').read_text())
    structure = PDBParser(QUIET=True).get_structure('3vhe', str(pdb))
    chain = structure[0]['A']
    ligand = chain[('H_42Q', 1170, ' ')]
    if contacts['structure_id'] != '3VHE' or contacts['ligand_comp_id'] != '42Q':
        raise BundleError('Expected the existing 3VHE / 42Q contact annotations')
    if len(list(ligand.get_atoms())) != contacts['n_ligand_heavy_atoms']:
        raise BundleError('Ligand atom count does not match the source annotations')
    mapping = {}
    for residue in chain:
        if residue.id[0] != ' ':
            continue
        key = f'3VHE:A:{residue.id[1]}{residue.id[2].strip()}'
        mapping[key] = {'chain': 'A', 'number': residue.id[1], 'insertion': residue.id[2].strip()}
    residues = []
    for contact in contacts['contacts']:
        # The upstream file omits insertion codes. Require a unique, blank-code
        # residue rather than guessing how to map an ambiguous author number.
        matches = [r for r in chain if r.id[0] == ' ' and r.id[1] == contact['residue_id']]
        if contact['chain'] != 'A' or len(matches) != 1 or matches[0].id[2] != ' ' or matches[0].resname != contact['residue_name']:
            raise BundleError('Ambiguous or mismatched source contact residue')
        residues.append({'key': f"3VHE:A:{contact['residue_id']}",
                         'label': f"{contact['residue_name']} A:{contact['residue_id']}",
                         'frequency': 1, 'core': True})
    unavailable = 'Structural reference only; no candidate scoring was performed.'
    manifest = {
        'schema_version': '1.0.0', 'run_id': '3vhe-reference', 'data_origin': 'reference_example',
        'disease': {'id': 'structural-reference', 'name': 'VEGFR2 binding example'},
        'target': {'id': 'P35968', 'name': 'VEGFR2 / KDR', 'structure_id': '3vhe-complex', 'site_kind': 'pocket',
                   'rationale': 'Experimental VEGFR2–42Q complex. This example shows a recorded binding pose, not a disease-treatment prediction.'},
        'signature': {'id': '3vhe-reference-site', 'source': 'known_ligand', 'status': 'available',
                      'reason': 'Existing 4.5 Å heavy-atom contact annotations for the bound 42Q ligand. Frequency 1 denotes presence in this single reference, not an ensemble.',
                      'threshold': 1, 'n_generated': 0, 'n_surviving': 0,
                      'addressability': {'status': 'not_evaluated', 'reason': 'No independent addressability assessment.'}, 'residues': residues},
        'structures': [{'id': '3vhe-complex', 'path': '3VHE_A_holo.pdb', 'sha256': hashlib.sha256(pdb.read_bytes()).hexdigest(),
                        'format': 'pdb', 'origin': 'experimental', 'frame_id': '3VHE-A-native',
                        'source': 'RCSB PDB 3VHE, chain A; X-ray diffraction, 1.55 Å. Local protein–ligand extract.',
                        'license': 'wwPDB archive data; cite PDB 3VHE and Oguro et al. (2010).',
                        'residue_map': mapping, 'ligand_residues': [{'chain': 'A', 'number': 1170, 'insertion': ''}]}],
        'candidates': [{'id': '42Q', 'name': '42Q · VEGFR2 inhibitor', 'modality': 'small_molecule', 'novelty': 'known_moa',
                        'rank': None, 'ranking_metric': 'Unranked experimental reference', 'structure_id': '3vhe-complex',
                        'approval': 'Research inhibitor; no regulatory approval is asserted.',
                        'coverage': {'status': 'not_evaluated', 'value': None, 'reason': unavailable}, 'engaged': [], 'missed': [],
                        'confidence': {'name': 'Prediction confidence', 'value': None, 'scale': 'Not applicable', 'gate': 'not_evaluated',
                                       'rule': 'Experimental structure; no prediction gate applied.', 'reason': unavailable},
                        'decoys': {'status': 'not_evaluated', 'reason': unavailable, 'percentile': None, 'scores': [], 'n': 0,
                                   'metric': 'core_coverage', 'context': 'P35968:3vhe-reference-site:small_molecule', 'tie_rule': 'Not applicable'},
                        'caveats': ['Crystal pose of the published inhibitor. Binding-site highlights reuse the pipeline’s geometric contact annotations; they are not a new prediction.'],
                        'source': 'https://www.rcsb.org/structure/3VHE'}],
        'evidence': [{'id': 'PDB:3VHE', 'title': 'VEGFR2 bound to inhibitor 42Q',
                      'claim': 'Experimental X-ray structure at 1.55 Å resolution.', 'url': 'https://www.rcsb.org/structure/3VHE'},
                     {'id': 'PMID:20833055', 'title': 'Oguro et al. (2010)',
                      'claim': 'Published synthesis and evaluation of pyrrolopyrimidine VEGFR2 inhibitors.',
                      'url': 'https://doi.org/10.1016/j.bmc.2010.08.017'}],
        'validation': [],
        'provenance': {'producer': 'visualization.reference', 'structure': 'RCSB PDB 3VHE', 'method': 'X-ray diffraction, 1.55 Å',
                       'binding_site': 'Upstream known_ligand_contacts.json; protein heavy atoms within 4.5 Å of ligand heavy atoms.',
                       'use': 'Reference display only. Excluded from prediction and hotspot-recovery evaluation.'},
    }
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copy2(pdb, destination / pdb.name)
    validate(manifest, destination)
    path = destination / 'manifest.json'
    path.write_text(json.dumps(manifest, indent=2))
    return path


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    print(make_reference(parser.parse_args().out))
