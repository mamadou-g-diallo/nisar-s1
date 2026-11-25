"""SNAP graph templates: every variable is substituted and the XML is well formed (no gpt needed)."""

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from nisar_s1.snap import render_graph

AOI = "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"
COMMON = dict(swath="IW2", aoi=AOI, region=AOI, out="/tmp/out.dim", pixel=20.0, crs="EPSG:32632")
VARS = {
    "s1_coherence": dict(ref="a.safe", sec="b.safe", cohWinRg=10, cohWinAz=3, **COMMON),
    "s1_coherence_radar": dict(ref="a.safe", sec="b.safe", pols="VV", cohWinRg=10, cohWinAz=3,
                               swath="IW2", aoi=AOI, region=AOI, out="/tmp/out.dim"),
    "s1_tc": dict(src="/tmp/radar.dim", pixel=20.0, crs="EPSG:32632", region=AOI, out="/tmp/out.dim"),
}


@pytest.mark.parametrize("name", sorted(VARS))
def test_graph_renders(tmp_path, name):
    xml = render_graph(name, tmp_path / f"{name}.xml", **VARS[name]).read_text()
    assert "${" not in xml
    root = ET.fromstring(xml)
    ops = [n.findtext("operator") for n in root.iter("node")]
    assert ops[-1] == "Write"


def test_two_step_split(tmp_path):
    """The radar step stops before terrain correction, the second step does only that."""
    radar = render_graph("s1_coherence_radar", Path(tmp_path / "r.xml"), **VARS["s1_coherence_radar"])
    ops = [n.findtext("operator") for n in ET.fromstring(radar.read_text()).iter("node")]
    assert "Coherence" in ops and "Terrain-Correction" not in ops and ops[-2] == "Subset"
    tc = render_graph("s1_tc", Path(tmp_path / "t.xml"), **VARS["s1_tc"])
    ops = [n.findtext("operator") for n in ET.fromstring(tc.read_text()).iter("node")]
    assert ops == ["Read", "Terrain-Correction", "Subset", "Write"]
