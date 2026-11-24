"""Georeferenced single-band arrays, GeoTIFF input/output and small helpers shared by all modules.

This module has no network dependency: it is the common currency between the NISAR reader
(:mod:`nisar_s1.io_nisar`), the Sentinel-1 reader (:mod:`nisar_s1.cband`), the Sentinel-2 reader
(:mod:`nisar_s1.optical`), land-cover masks (:mod:`nisar_s1.landcover`) and the analysis modules.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform_bounds

# GDAL settings for fast remote access to cloud-optimised GeoTIFFs (fewer HTTP requests).
GDAL_ENV = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MULTIRANGE="YES",
                GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES", VSI_CACHE="TRUE",
                CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff,.TIF", GDAL_CACHEMAX=512)


@dataclass
class Raster:
    """A single-band array with its georeferencing (``transform`` maps pixel corners)."""
    data: np.ndarray
    transform: rasterio.Affine
    crs: rasterio.crs.CRS

    @property
    def pixel_area_m2(self) -> float:
        return abs(self.transform.a * self.transform.e)

    @property
    def res(self) -> float:
        """Pixel size along x (projection units, metres for UTM)."""
        return abs(self.transform.a)

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        """(xmin, ymin, xmax, ymax) of the pixel edges."""
        h, w = self.data.shape
        x0, y0 = self.transform * (0, 0)
        x1, y1 = self.transform * (w, h)
        return min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)


def bbox_to_crs(bbox_wgs84, crs) -> tuple[float, float, float, float]:
    """Project a ``[lon_min, lat_min, lon_max, lat_max]`` box to ``crs`` (densified edges, so the
    result encloses the whole box)."""
    return tuple(transform_bounds("EPSG:4326", crs, *bbox_wgs84, densify_pts=21))


def to_db(linear: np.ndarray) -> np.ndarray:
    """10*log10 of a power array; zero, negative or missing values become NaN."""
    a = np.asarray(linear, dtype="float64")
    with np.errstate(divide="ignore", invalid="ignore"):
        return (10 * np.log10(np.where(a > 0, a, np.nan))).astype("float32")


def save_raster(r: Raster, path, description: str | None = None) -> None:
    """Write a :class:`Raster` as a compressed float32 GeoTIFF (NaN = no data)."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", driver="GTiff", height=r.data.shape[0], width=r.data.shape[1],
                       count=1, dtype="float32", crs=r.crs, transform=r.transform, nodata=np.nan,
                       compress="deflate", predictor=3, tiled=True) as dst:
        dst.write(r.data.astype("float32"), 1)
        if description:
            dst.set_band_description(1, description)


def load_raster(path) -> Raster:
    """Read band 1 of a GeoTIFF into a :class:`Raster`."""
    with rasterio.open(path) as src:
        return Raster(src.read(1), src.transform, src.crs)
