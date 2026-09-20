import base64
import html
import json
from pathlib import Path


def render_report(bundle, renders=None):
    escape = lambda value: html.escape(str(value), quote=True)
    m = bundle["manifest"]
    origin = m["data_origin"].replace("_", " ").upper()
    parts = ["<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width'><title>Rebind — report</title><style>body{font:16px/1.6 system-ui,sans-serif;color:#172b3a;max-width:1080px;margin:48px auto;padding:0 28px}h1{font-size:42px;line-height:1.15}h2{margin-top:48px;color:#007f73}small{color:#52616b}table{border-collapse:collapse;width:100%;font-size:14px}td,th{text-align:left;padding:12px;border-bottom:1px solid #dbe4e8}aside{padding:16px;border:1px solid #b45309;background:#fff8eb}.metric{font-size:28px}.muted{color:#64748b}img{max-width:100%}.badge{font-weight:700;letter-spacing:.08em}code{overflow-wrap:anywhere}@media print{body{margin:0}h2{break-after:avoid}tr,figure{break-inside:avoid}}</style>"]
    parts += [f"<p class='badge'>REBIND / {escape(m['run_id'])}</p><h1>{escape(m['disease']['name'])}<br><small>{escape(m['target']['name'])}</small></h1>",
              f"<aside><b>{escape(origin)}</b><br>{escape(bundle['disclaimer'])}</aside>",
              f"<h2>01 / Evidence & target</h2><p>{escape(m['target']['rationale'])}</p>"]
    for e in m["evidence"]:
        parts.append(f"<h3>{escape(e['title'])}</h3><p>{escape(e['claim'])}</p><small>{escape(e['id'])}</small>")
        if e["url"]:
            parts.append(f" <a href='{escape(e['url'])}' rel='noopener noreferrer'>Source</a>")
    sig = m["signature"]
    parts.append(f"<h2>02 / Interface signature</h2><p>Source: <b>{escape(sig['source'])}</b>. Status: {escape(sig['status'])}. {escape(sig['reason'])}</p><p>Ensemble counts: {sig['n_surviving']} surviving / {sig['n_generated']} generated. These counts are not applicable to a geometry-only or known-ligand signature unless explicitly populated upstream.</p><p>Small-molecule addressability: <b>{escape(sig['addressability']['status'])}</b>. {escape(sig['addressability']['reason'])}</p>")
    parts.append("<table><thead><tr><th>Residue</th><th>Engagement frequency</th><th>Core</th></tr></thead><tbody>")
    for r in sig["residues"]:
        parts.append(f"<tr><td>{escape(r['label'])}</td><td>{r['frequency']:.0%}</td><td>{'Yes' if r['core'] else 'No'}</td></tr>")
    parts.append("</tbody></table><h2>03 / Candidate hypotheses</h2>")
    for modality in ("small_molecule", "peptide", "biologic"):
        candidates = [c for c in m["candidates"] if c["modality"] == modality]
        parts.append(f"<h3>{escape(modality.replace('_', ' ').title())}</h3>")
        suppressed = modality == "small_molecule" and sig["addressability"]["status"] != "addressable"
        if suppressed:
            parts.append(f"<aside>Leaderboard suppressed: {escape(sig['addressability']['reason'])}. Records below are inspection-only.</aside>")
        parts.append(f"<p>n = {len(candidates)} records; {sum(c['rank'] is not None for c in candidates)} upstream-ranked.</p>")
        parts.append("<table><thead><tr><th>Rank / candidate</th><th>Coverage</th><th>Confidence gate</th><th>Calibration</th><th>Missed core / caveats</th></tr></thead><tbody>")
        for c in sorted(candidates, key=lambda c: (c["rank"] is None, c["rank"] or 0)):
            coverage = f"{c['coverage']['value']:.0%}" if c["coverage"]["value"] is not None else "Not evaluated"
            d = c["decoys"]
            calibrated = f"{d['percentile']:g} percentile; n={d['n']}" if d["status"] == "available" else "Uncalibrated: " + d["reason"]
            parts.append(f"<tr><td>{escape(c['rank'] if not suppressed and c['rank'] is not None else '—')} / {escape(c['name'])}<br><small>{escape(c['novelty'])}; {escape(c['approval'])}</small></td><td>{coverage}</td><td>{escape(c['confidence']['gate'])}<br><small>{escape(c['confidence']['name'])}: {escape(c['confidence']['value'])}; {escape(c['confidence']['scale'])}</small></td><td>{escape(calibrated)}</td><td>{escape(', '.join(c['missed']) or 'None reported')}<br>{escape('; '.join(c['caveats']))}</td></tr>")
        parts.append("</tbody></table>")
    parts.append("<h2>04 / Validation & limitations</h2>")
    for v in m["validation"]:
        parts.append(f"<p><b>{escape(v['label'])}</b>: {escape(v['status'])}; n={v['n']}; value={escape(v['value'] if v['value'] is not None else 'not evaluated')}. {escape(v['reason'])}<br><small>Source: {escape(v['source'])}</small></p>")
    parts.append("<h2>05 / Molecular figures</h2>")
    if renders:
        from .bundle import BundleError, safe_asset
        renders = Path(renders)
        record = json.loads((renders / "render_manifest.json").read_text())
        if record["run_id"] != m["run_id"] or record["origin"] != m["data_origin"]:
            raise BundleError("Render manifest does not match run/origin")
        for scene in record["scenes"]:
            path = safe_asset(renders, scene["file"])
            raw = path.read_bytes()
            if not raw.startswith(b"\x89PNG\r\n\x1a\n"):
                raise BundleError("Expected rendered PNG")
            data = base64.b64encode(raw).decode()
            parts.append(f"<figure><img src='data:image/png;base64,{data}' alt='{escape(scene['caption'])}'><figcaption>{escape(origin)} — {escape(scene['caption'])}</figcaption></figure>")
    else:
        parts.append("<p>Molecular stills not rendered. Use the interactive bundle for structures; a verified PyMOL installation is required to generate presentation figures.</p>")
    parts.append("<h2>06 / Provenance</h2><table>")
    for key, value in m["provenance"].items():
        parts.append(f"<tr><th>{escape(key)}</th><td>{escape(value)}</td></tr>")
    parts.append("</table>")
    for warning in bundle["warnings"]:
        parts.append(f"<p>{escape(warning)}</p>")
    parts.append(f"<p><small>{escape(bundle['disclaimer'])} This snapshot makes no significance or superiority claim without the corresponding scientific analysis.</small></p></html>")
    return "".join(parts)
