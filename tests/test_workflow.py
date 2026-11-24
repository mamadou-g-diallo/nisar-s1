"""Offline tests of the generic workflow on synthetic NISAR products shaped like the real ones
in UTM zone 28N (40 MHz mode: 10 m GCOV + 80 m side band B, dual- or single-pol,
GUNW coherence at 20 m)."""

import h5py
import numpy as np
import pandas as pd
import pytest
import rasterio
from affine import Affine
from pyproj import Transformer

from nisar_s1 import cli, io_nisar
from nisar_s1.cband import read_cog
from nisar_s1.landcover import class_separability, class_stats, erode_classes
from nisar_s1.raster import Raster, load_raster, to_db
from nisar_s1.workflow import catalog_summary, extract_nisar, summarize

EPSG = 32628
X0, Y0 = Transformer.from_crs(4326, EPSG, always_xy=True).transform(-16.62, 28.32)
X0, Y0 = np.floor(X0 / 1000) * 1000 - 1000, np.ceil(Y0 / 1000) * 1000 + 1000  # grid corner
BBOX = [-16.62, 28.28, -16.58, 28.32]  # about 4.3 x 4.4 km, inside the products


def _grid_group(g, n, step):
    g["xCoordinates"] = X0 + step / 2 + step * np.arange(n)
    g["yCoordinates"] = Y0 - step / 2 - step * np.arange(n)
    g["projection"] = np.uint32(EPSG)


def make_gcov(path, pols=("HH", "HV")):
    """10 m frequency A (checkerboard 0.1 / 0.3 in HH) and 80 m frequency B."""
    with h5py.File(path, "w") as h:
        h["science/LSAR/identification/listOfFrequencies"] = np.array([b"A", b"B"])
        for f, step, n in (("A", 10.0, 800), ("B", 80.0, 100)):
            g = h.create_group(f"science/LSAR/GCOV/grids/frequency{f}")
            _grid_group(g, n, step)
            g["listOfPolarizations"] = np.array([p.encode() for p in pols])
            ii, jj = np.indices((n, n))
            hh = np.where((ii + jj) % 2 == 0, 0.1, 0.3).astype("float32")
            g["HHHH"] = hh
            if "HV" in pols:
                g["HVHV"] = np.full((n, n), 0.01, "float32")
            g["rtcGammaToSigmaFactor"] = np.full((n, n), 1.0, "float32")


def make_gunw(path, coh=0.6):
    with h5py.File(path, "w") as h:
        w = h.create_group("science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram")
        _grid_group(w, 400, 20.0)
        h["science/LSAR/GUNW/grids/frequencyA/listOfPolarizations"] = np.array([b"HH"])
        w["HH/coherenceMagnitude"] = np.full((400, 400), coh, "float32")


def test_polarizations_frequencies_and_bbox(tmp_path):
    make_gcov(tmp_path / "dual.h5")
    make_gcov(tmp_path / "single.h5", pols=("HH",))
    with io_nisar.open_product(tmp_path / "dual.h5") as h:
        assert io_nisar.frequencies(h) == ["A", "B"]
        assert io_nisar.polarizations(h, "GCOV") == ["HH", "HV"]
        assert io_nisar.posting(h, "GCOV", frequency="B") == 80
        b = io_nisar.gcov_backscatter(h, bbox_wgs84=BBOX)
        assert set(b) == {"HH", "HV"} and b["HH"].res == 10
        assert np.allclose(b["HV"].data, -20)
        bb = io_nisar.gcov_backscatter(h, bbox_wgs84=BBOX, frequency="B")
        assert bb["HH"].res == 80
    with io_nisar.open_product(tmp_path / "single.h5") as h:
        assert set(io_nisar.gcov_backscatter(h, bbox_wgs84=BBOX)) == {"HH"}
        with pytest.raises(KeyError):
            io_nisar.gcov_backscatter(h, bbox_wgs84=BBOX, pols=["HV"])


def _catalog(tmp_path):
    rows = []
    for k, d in enumerate(("2026-09-01", "2026-09-13")):
        p = tmp_path / f"gcov{k}.h5"
        make_gcov(p)
        rows.append({"product": "GCOV", "track": 124, "direction": "D", "frame": 82,
                     "bw_code": "4005", "pol_code": "DHDH", "sceneName": p.name, "url": str(p),
                     "ref_date": pd.Timestamp(d).date(), "startTime": pd.Timestamp(d, tz="UTC")})
    p = tmp_path / "gunw.h5"
    make_gunw(p)
    rows.append({"product": "GUNW", "track": 124, "direction": "D", "frame": 82, "bw_code": "4000",
                 "pol_code": "SH", "sceneName": p.name, "url": str(p),
                 "ref_date": pd.Timestamp("2026-09-01").date(),
                 "sec_date": pd.Timestamp("2026-09-13").date(), "baseline_days": 12,
                 "startTime": pd.Timestamp("2026-09-01", tz="UTC")})
    return pd.DataFrame(rows)


def test_extract_on_common_grid_and_resume(tmp_path):
    cat = _catalog(tmp_path)
    out = tmp_path / "out"
    idx = extract_nisar(BBOX, "2026-09-01", "2026-09-30", out, "GCOV", catalog=cat, say=lambda m: None)
    idx = extract_nisar(BBOX, "2026-09-01", "2026-09-30", out, "GUNW", catalog=cat, say=lambda m: None)
    assert set(idx["variable"]) == {"L HH (dB)", "L HV (dB)", "L coh20"}
    rasters = [load_raster(p) for p in idx["path"]]
    assert len({r.data.shape for r in rasters}) == 1 and rasters[0].res == 20
    hh = load_raster(idx[idx.variable == "L HH (dB)"]["path"].iloc[0]).data
    # 2x2 average of the 0.1 / 0.3 checkerboard in linear power = 0.2 (not the mean of the dB)
    assert np.nanmedian(hh) == pytest.approx(10 * np.log10(0.2), abs=0.01)
    coh = load_raster(idx[idx.variable == "L coh20"]["path"].iloc[0]).data
    assert np.nanmedian(coh) == pytest.approx(0.6, abs=1e-4)
    assert len(catalog_summary(cat)) == 2

    def no_network(url):
        raise AssertionError("already extracted products must not be read again")
    idx2 = extract_nisar(BBOX, "2026-09-01", "2026-09-30", out, "GCOV", catalog=cat,
                         open_fn=no_network, say=lambda m: None)
    assert len(idx2) == len(idx)


def test_summarize_with_given_classes(tmp_path):
    cat = _catalog(tmp_path)
    out = tmp_path / "out"
    idx = extract_nisar(BBOX, "2026-09-01", "2026-09-30", out, "GCOV", catalog=cat, say=lambda m: None)
    shape = load_raster(idx["path"].iloc[0]).data.shape
    classes = np.full(shape, 40, "uint8")
    classes[:, : shape[1] // 2] = 10
    stats, sep = summarize(out, classes=classes, say=lambda m: None)
    assert set(stats["class"]) == {10, 40}
    assert (out / "class_profiles.png").exists() and (out / "separability.csv").exists()
    assert set(sep["class_b"]) == {40}


def test_class_tools():
    tr = Affine(20, 0, 0, 0, -20, 0)
    classes = np.zeros((40, 40), "uint8")
    classes[:, :20], classes[:, 20:] = 10, 95
    a = np.where(classes == 10, -12.0, -18.0).astype("float32")
    a += np.random.default_rng(0).normal(0, 0.5, a.shape).astype("float32")
    layers = {"L HV (dB)": Raster(a, tr, "EPSG:32628")}
    st = class_stats(layers, classes)
    assert st.set_index("class").loc[10, "median"] == pytest.approx(-12, abs=0.2)
    eroded = erode_classes(classes, 1)
    assert (eroded[:, 19] == 0).all() and (eroded[:, 20] == 0).all() and (eroded[:, 10] == 10).all()
    sep = class_separability(layers, classes, 10, 95)
    assert sep["auc_sep"].iloc[0] > 0.99


def test_read_cog_average_in_linear(tmp_path):
    p = tmp_path / "s1.tif"
    a = np.where(np.indices((40, 40)).sum(0) % 2 == 0, 0.1, 0.3).astype("float32")
    a[0, 0] = -32768
    with rasterio.open(p, "w", driver="GTiff", height=40, width=40, count=1, dtype="float32",
                       crs="EPSG:32628", transform=Affine(10, 0, 300000, 0, -10, 1500000),
                       nodata=-32768) as ds:
        ds.write(a, 1)
    grid = Raster(np.full((20, 20), np.nan, "float32"), Affine(20, 0, 300000, 0, -20, 1500000),
                  rasterio.crs.CRS.from_epsg(32628))
    r = read_cog(str(p), grid=grid)
    assert r.data.shape == (20, 20)
    assert np.nanmedian(r.data) == pytest.approx(0.2, abs=1e-6)
    assert np.isfinite(to_db(r.data)).sum() >= 399


def test_cli_parsing(capsys):
    ap = cli.build_parser()
    a = ap.parse_args(["search", "--bbox", "6.70", "47.85", "6.82", "47.95",
                       "--start", "2026-07-01", "--end", "2026-11-30"])
    assert cli._bbox(a) == [6.70, 47.85, 6.82, 47.95]
    a = ap.parse_args(["extract", "--site", "vosges", "--start", "2026-07-01", "--end",
                       "2026-11-30", "--out", "x", "--product", "GCOV"])
    assert cli._bbox(a) == cli.EXAMPLE_SITES["vosges"][0] and a.product == ["GCOV"]
    cli.main(["sites"])
    assert "vosges" in capsys.readouterr().out


def test_user_sites_file(tmp_path, monkeypatch, capsys):
    f = tmp_path / "sites.yaml"
    f.write_text("demo:\n  bbox: [10.0, 45.0, 10.2, 45.2]\n  label: my test site\n")
    ap = cli.build_parser()
    a = ap.parse_args(["search", "--sites-file", str(f), "--site", "demo",
                       "--start", "2026-07-01", "--end", "2026-11-30"])
    assert cli._bbox(a) == [10.0, 45.0, 10.2, 45.2]
    monkeypatch.setenv("NISAR_S1_SITES", str(f))
    cli.main(["sites"])
    out = capsys.readouterr().out
    assert "demo" in out and "my test site" in out and "vosges" in out
    f.write_text("bad:\n  bbox: [1, 2, 3]\n")
    with pytest.raises(SystemExit):
        cli.load_sites(f)


def test_shared_geometry_prefers_a_track_with_both_products():
    from nisar_s1.workflow import shared_geometry
    rows = [("GCOV", 131, "A", "DHDH")] * 5 + [("GCOV", 23, "D", "DHDH")] * 3 + \
           [("GUNW", 23, "D", "SH")] * 4
    df = pd.DataFrame(rows, columns=["product", "track", "direction", "pol_code"])
    assert shared_geometry(df) == (23, "D")
    assert shared_geometry(df[df["product"] == "GCOV"]) == (131, "A")


def test_index_independent_of_working_directory(tmp_path, monkeypatch):
    """An output folder given as a relative path stays readable from anywhere and after a move."""
    from nisar_s1.workflow import read_index, temporal_medians
    cat = _catalog(tmp_path)
    monkeypatch.chdir(tmp_path)
    extract_nisar(BBOX, "2026-09-01", "2026-09-30", "out", "GCOV", catalog=cat, say=lambda m: None)
    assert not any("/" in p for p in pd.read_csv(tmp_path / "out" / "files.csv")["path"])
    (tmp_path / "out").rename(tmp_path / "moved")
    monkeypatch.chdir("/")
    idx = read_index(tmp_path / "moved")
    assert all(load_raster(p).data.size for p in idx["path"])
    assert len(temporal_medians(tmp_path / "moved")) == 2
