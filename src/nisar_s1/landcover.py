"""Land-cover reference (ESA WorldCover 2021, 10 m, public COGs) and per-class statistics.

The typical first look at a new site is: which land-cover classes does each SAR variable separate,
and how differently in L band (NISAR) and C band (Sentinel-1)? :func:`worldcover_like` puts the
WorldCover classes on the grid of any raster, :func:`class_stats` summarises each variable per
class and :func:`class_separability` gives the separation AUC between two classes.

WorldCover is a map, not ground truth (overall accuracy about 75 % in 2021, Zanaga et al. 2022):
class statistics are indicative, and boundaries should be eroded before comparing classes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling

from .raster import GDAL_ENV, Raster
from .stats import separability

WORLDCOVER_URL = ("https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
                  "ESA_WorldCover_10m_2021_v200_{tile}_Map.tif")
WORLDCOVER_TREE = 10
WORLDCOVER_CLASSES = {
    10: "Tree cover", 20: "Shrubland", 30: "Grassland", 40: "Cropland", 50: "Built-up",
    60: "Bare / sparse vegetation", 70: "Snow and ice", 80: "Permanent water bodies",
    90: "Herbaceous wetland", 95: "Mangroves", 100: "Moss and lichen",
}


def _worldcover_tiles(bounds_wgs84) -> list[str]:
    """Names of the 3x3 degree WorldCover tiles intersecting ``bounds_wgs84``."""
    x0, y0, x1, y1 = bounds_wgs84
    names = []
    for lat in range(int(np.floor(y0 / 3) * 3), int(np.floor(y1 / 3) * 3) + 1, 3):
        for lon in range(int(np.floor(x0 / 3) * 3), int(np.floor(x1 / 3) * 3) + 1, 3):
            names.append(f"{'N' if lat >= 0 else 'S'}{abs(lat):02d}"
                         f"{'E' if lon >= 0 else 'W'}{abs(lon):03d}")
    return names


def worldcover_like(ref: Raster, bounds_wgs84=None) -> np.ndarray:
    """WorldCover class codes (uint8, 0 = no data) on the grid of ``ref`` (mode resampling).

    ``bounds_wgs84`` only selects the tiles to read; by default it is derived from ``ref``.
    """
    from rasterio.vrt import WarpedVRT
    from rasterio.warp import transform_bounds
    if bounds_wgs84 is None:
        bounds_wgs84 = transform_bounds(ref.crs, "EPSG:4326", *ref.bounds, densify_pts=21)
    out = np.zeros(ref.data.shape, dtype="uint8")
    h, w = ref.data.shape
    for tile in _worldcover_tiles(bounds_wgs84):
        with rasterio.Env(**GDAL_ENV), rasterio.open(WORLDCOVER_URL.format(tile=tile)) as src, \
                WarpedVRT(src, crs=ref.crs, transform=ref.transform, width=w, height=h,
                          resampling=Resampling.mode, nodata=0) as vrt:
            a = vrt.read(1)
            out = np.where(a > 0, a, out)
    return out


def erode_classes(classes: np.ndarray, pixels: int = 1) -> np.ndarray:
    """Set to 0 every pixel within ``pixels`` of a different class (removes mixed edge pixels)."""
    from scipy import ndimage
    out = classes.copy()
    for c in np.unique(classes):
        if c == 0:
            continue
        m = classes == c
        core = ndimage.binary_erosion(m, iterations=pixels, border_value=1)
        out[m & ~core] = 0
    return out


def class_stats(layers: dict[str, Raster], classes: np.ndarray, names: dict | None = None,
                min_pixels: int = 50) -> pd.DataFrame:
    """Median, quartiles, mean and pixel count of every layer in every class.

    ``layers`` maps a variable name (e.g. ``"L HV (dB)"``) to a :class:`Raster` on the same grid
    as ``classes``. Classes with fewer than ``min_pixels`` valid pixels are dropped.
    """
    names = WORLDCOVER_CLASSES if names is None else names
    rows = []
    for var, r in layers.items():
        if r.data.shape != classes.shape:
            raise ValueError(f"{var}: shape {r.data.shape} differs from the class map {classes.shape}")
        for c in np.unique(classes):
            if c == 0:
                continue
            v = r.data[(classes == c) & np.isfinite(r.data)]
            if v.size < min_pixels:
                continue
            q25, q50, q75 = np.percentile(v, [25, 50, 75])
            rows.append({"variable": var, "class": int(c), "class_name": names.get(int(c), str(c)),
                         "n": int(v.size), "median": float(q50), "q25": float(q25),
                         "q75": float(q75), "mean": float(v.mean())})
    return pd.DataFrame(rows)


def class_separability(layers: dict[str, Raster], classes: np.ndarray, class_a: int, class_b: int,
                       max_samples: int = 20000, seed: int = 0) -> pd.DataFrame:
    """Separation of ``class_a`` from ``class_b`` for every layer (AUC, Cohen's d, ...).

    At most ``max_samples`` pixels per class are drawn at random (fixed seed), which keeps the
    Mann-Whitney AUC fast; neighbouring pixels are correlated, so the AUC is a descriptive effect
    size, not a test.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for var, r in layers.items():
        a = r.data[(classes == class_a) & np.isfinite(r.data)]
        b = r.data[(classes == class_b) & np.isfinite(r.data)]
        a = rng.choice(a, max_samples, replace=False) if a.size > max_samples else a
        b = rng.choice(b, max_samples, replace=False) if b.size > max_samples else b
        rows.append({"variable": var, "class_a": class_a, "class_b": class_b, **separability(a, b)})
    return pd.DataFrame(rows)
