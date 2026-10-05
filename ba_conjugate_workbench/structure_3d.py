"""3Dmol/py3Dmol rendering helpers."""

from __future__ import annotations

import html
import json

from .highlighting import infer_highlights


def molblock_to_html(mol_block: str, row: dict, label_mode: str = "Off", height: int = 520) -> str:
    """Return an embeddable py3Dmol HTML viewer.

    3Dmol.js provides left-drag rotation, right-drag translation/pan, and
    wheel zoom in the browser.
    """
    safe_block = mol_block.replace("\\", "\\\\").replace("`", "\\`")
    spec = infer_highlights(row)
    label_script = ""
    if label_mode == "Element only":
        label_script = "viewer.addPropertyLabels('elem', {}, {fontSize: 10, fontColor: 'black', showBackground: false});"
    elif label_mode == "Element + RDKit atom index":
        label_script = "viewer.addPropertyLabels('index', {}, {fontSize: 10, fontColor: 'black', showBackground: false});"
    elif label_mode == "BA conventional labels + atom indices":
        label_map = json.dumps({str(k): f"{v} (idx {k})" for k, v in spec.labels.items()})
        label_script = f"""
const baLabels = {label_map};
viewer.getModel().selectedAtoms({{}}).forEach(function(atom) {{
  const text = baLabels[String(atom.index)];
  if (text) {{
    viewer.addLabel(text, {{
      position: {{x: atom.x, y: atom.y, z: atom.z}},
      fontSize: 11,
      fontColor: "black",
      backgroundColor: "white",
      backgroundOpacity: 0.75
    }});
  }}
}});
"""
    site_atoms = ",".join(str(i) for i in spec.site_atoms)
    partner_atoms = ",".join(str(i) for i in spec.partner_atoms[:200])
    core_atoms = ",".join(str(i) for i in spec.core_atoms[:200])
    return f"""
<div id="ba_viewer" style="width: 100%; height: {height}px; position: relative;"></div>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
const element = document.getElementById("ba_viewer");
const viewer = $3Dmol.createViewer(element, {{backgroundColor: "white"}});
viewer.addModel(`{safe_block}`, "sdf");
viewer.setStyle({{}}, {{stick: {{radius: 0.16, colorscheme: "Jmol"}}}});
if ("{core_atoms}".length) {{
  viewer.setStyle({{index: [{core_atoms}]}}, {{stick: {{radius: 0.18, color: "#4f7fd8"}}}});
}}
if ("{partner_atoms}".length) {{
  viewer.setStyle({{index: [{partner_atoms}]}}, {{stick: {{radius: 0.18, color: "#dd8a23"}}}});
}}
if ("{site_atoms}".length) {{
  viewer.setStyle({{index: [{site_atoms}]}}, {{stick: {{radius: 0.26, color: "#df2f2f"}}, sphere: {{scale: 0.32, color: "#df2f2f"}}}});
}}
{label_script}
viewer.zoomTo();
viewer.render();
</script>
"""


def escaped_molblock(mol_block: str) -> str:
    return html.escape(mol_block)
