"""Sentinel-1 RTC helpers that need no network: orbit choice and per-pixel composite."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import rasterio
from affine import Affine

from nisar_s1 import cband
from nisar_s1.grid import make_grid


def test_common_orbit_most_periods():
    by_year = {2019: pd.DataFrame({"relative_orbit": [91, 91, 18]}),
               2021: pd.DataFrame({"relative_orbit": [91]}),
               2023: pd.DataFrame({"relative_orbit": [18, 91]}),
               2025: pd.DataFrame({"relative_orbit": []})}
    assert cband.common_orbit(by_year) == 91


def _cog(path, values, nodata=-32768.0):
    a = np.asarray(values, "float32")
    with rasterio.open(path, "w", driver="GTiff", height=a.shape[0], width=a.shape[1], count=1,
                       dtype="float32", crs="EPSG:32619", transform=Affine(10, 0, 300000, 0, -10, 5000040),
                       nodata=nodata) as ds:
        ds.write(a, 1)
    return SimpleNamespace(assets={"vv": SimpleNamespace(href=str(path))})


def test_rtc_composite_median_ignores_nodata(tmp_path, monkeypatch):
    monkeypatch.setattr(cband, "signed_href", lambda href, collection=cband.S1_RTC: href)
    nd = -32768.0
    items = [_cog(tmp_path / "a.tif", [[1, 1, 1, 1]] * 4),
             _cog(tmp_path / "b.tif", [[3, 3, 3, 3]] * 4),
             _cog(tmp_path / "c.tif", [[2, 2, nd, nd]] * 4)]
    grid = make_grid((300000, 5000000, 300040, 5000040), 20, "EPSG:32619", snap=False)
    out = cband.rtc_composite(items, "vv", grid, resampling="nearest")
    assert out.data.shape == (2, 2)
    np.testing.assert_allclose(out.data[:, 0], 2.0)   # median of 1, 3, 2
    np.testing.assert_allclose(out.data[:, 1], 2.0)   # median of 1, 3 (third scene is no-data)


def test_make_grid_without_snapping_starts_at_corner():
    g = make_grid((300005, 5000003, 300046, 5000043), 20, "EPSG:32619", snap=False)
    assert (g.transform.c, g.transform.f) == (300005, 5000043)
    assert g.data.shape == (2, 3)                     # ceil(40/20), ceil(41/20)
    s = make_grid((300005, 5000003, 300046, 5000043), 20, "EPSG:32619")
    assert (s.transform.c, s.transform.f) == (300000, 5000060)
