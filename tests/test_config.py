"""Configuration lookup: explicit path, NISAR_S1_ROOT, missing files (no network)."""
from pathlib import Path

import pytest

from nisar_s1 import config

SITE_YAML = """nom: test
chemins:
  results: results/test
  data: ~/donnees_test
"""
PARAMS_YAML = """aoi_recherche:
  bbox_wgs84: [6.7, 47.8, 6.9, 48.0]
chemins:
  data: ~/donnees
"""


def _make_root(tmp_path: Path) -> Path:
    (tmp_path / "config" / "sites").mkdir(parents=True)
    (tmp_path / "config" / "sites" / "test.yaml").write_text(SITE_YAML)
    (tmp_path / "config" / "params.yaml").write_text(PARAMS_YAML)
    return tmp_path


def test_load_site_from_explicit_yaml_path(tmp_path):
    root = _make_root(tmp_path)
    site = config.load_site(root / "config" / "sites" / "test.yaml", root=root)
    assert site["chemins"]["results"] == str(root / "results" / "test")
    assert site["chemins"]["data"] == str(Path("~/donnees_test").expanduser())


def test_root_from_environment_variable(tmp_path, monkeypatch):
    root = _make_root(tmp_path)
    monkeypatch.setenv("NISAR_S1_ROOT", str(root))
    assert config.project_root() == root.resolve()
    assert config.load_site("test")["chemins"]["results"] == str(root.resolve() / "results" / "test")
    params = config.load_params()
    assert config.search_aoi_wkt(params).startswith("POLYGON((6.7 47.8")


def test_missing_file_gives_a_clear_error(tmp_path, monkeypatch):
    monkeypatch.setenv("NISAR_S1_ROOT", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="NISAR_S1_ROOT"):
        config.load_site("absent")


def test_bbox_to_wkt_is_closed():
    wkt = config.bbox_to_wkt([0, 1, 2, 3])
    assert wkt == "POLYGON((0 1,2 1,2 3,0 3,0 1))"


def test_load_kit_finds_training_layout(tmp_path):
    from nisar_s1.config import load_kit
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "params.yaml").write_text(f"chemins:\n  data_root: {tmp_path / 'data'}\n")
    nb = tmp_path / "training" / "notebooks"
    nb.mkdir(parents=True)
    (nb / "kit_config.yaml").write_text(
        "data_dir: null\ngcov:\n  granule: G1\n  url: https://x/{granule}.h5\n"
        "aoi:\n  site_polygons: config/site.geojson\n")
    kit = load_kit(root=tmp_path)
    assert kit["gcov"]["url"] == "https://x/G1.h5"
    assert kit["aoi"]["site_polygons"] == str(tmp_path / "config" / "site.geojson")
