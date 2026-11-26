"""C-band reference without SNAP: Sentinel-1 radiometrically terrain-corrected backscatter.

Source: the ``sentinel-1-rtc`` collection of Microsoft Planetary Computer (global, IW GRD
processed to terrain-flattened gamma0, linear power, 10 m, UTM cloud-optimised GeoTIFFs). Access
is anonymous: a short-lived read token is requested from the public token endpoint, so no account
is needed. This gives the C-band backscatter counterpart of NISAR GCOV on any site in a few lines;
Sentinel-1 *coherence* still requires SLC processing (see :mod:`nisar_s1.snap`).
"""

from __future__ import annotations

import time
import warnings

import numpy as np
import pandas as pd
import rasterio
import requests
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import from_bounds

from .raster import GDAL_ENV, Raster, bbox_to_crs, to_db

PC_STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
PC_TOKEN = "https://planetarycomputer.microsoft.com/api/sas/v1/token/{collection}"
S1_RTC = "sentinel-1-rtc"
_NODATA = -32768.0
_tokens: dict[str, tuple[str, float]] = {}


def _token(collection: str = S1_RTC) -> str:
    """Anonymous read token for a Planetary Computer collection (cached until 5 min before expiry)."""
    tok, expiry = _tokens.get(collection, ("", 0.0))
    if time.time() < expiry - 300:
        return tok
    r = requests.get(PC_TOKEN.format(collection=collection), timeout=30)
    r.raise_for_status()
    j = r.json()
    exp = pd.Timestamp(j["msft:expiry"]).timestamp()
    _tokens[collection] = (j["token"], exp)
    return j["token"]


def signed_href(href: str, collection: str = S1_RTC) -> str:
    """Append the read token to an asset URL."""
    return f"{href}{'&' if '?' in href else '?'}{_token(collection)}"


def search_s1_rtc(bbox_wgs84, start: str, end: str, orbit_state: str | None = None,
                  relative_orbit: int | None = None, max_items: int = 500) -> pd.DataFrame:
    """Sentinel-1 RTC scenes intersecting ``bbox_wgs84`` between ``start`` and ``end``.

    One row per scene: id, datetime, platform, orbit_state, relative_orbit, item (the STAC item).
    """
    from pystac_client import Client
    items = Client.open(PC_STAC).search(collections=[S1_RTC], bbox=list(bbox_wgs84),
                                        datetime=f"{start}/{end}", max_items=max_items).items()
    rows = [{"id": i.id, "datetime": pd.Timestamp(i.properties["datetime"]),
             "platform": i.properties.get("platform"),
             "orbit_state": i.properties.get("sat:orbit_state"),
             "relative_orbit": i.properties.get("sat:relative_orbit"), "item": i} for i in items]
    df = pd.DataFrame(rows)
    if not len(df):
        return df
    if orbit_state:
        df = df[df.orbit_state == orbit_state]
    if relative_orbit is not None:
        df = df[df.relative_orbit == relative_orbit]
    return df.sort_values("datetime").reset_index(drop=True)


def common_orbit(scenes_by_period: dict) -> int:
    """Relative orbit present in the largest number of periods (years, seasons...).

    ``scenes_by_period`` maps a period to the table returned by :func:`search_s1_rtc`. Using one
    relative orbit for every period keeps the incidence angle of each pixel constant, so that a
    change between periods is not a change of viewing geometry."""
    counts: dict[int, int] = {}
    for df in scenes_by_period.values():
        for o in set(df.relative_orbit) if len(df) else ():
            counts[o] = counts.get(o, 0) + 1
    if not counts:
        raise ValueError("no scene in any period")
    return max(counts, key=counts.get)   # first orbit with the highest count


def rtc_composite(items, pol: str, grid: Raster, stat: str = "median",
                  resampling: str = "bilinear") -> Raster:
    """Per-pixel composite (``median`` or ``mean``) of several Sentinel-1 RTC scenes on ``grid``.

    Each scene's ``pol`` asset (``"vv"`` or ``"vh"``, linear gamma0) is reprojected onto the grid,
    then combined ignoring no-data. The result stays in linear power (NaN = no data). A median over
    the scenes of a season reduces speckle and short-lived moisture effects."""
    from rasterio.warp import reproject

    h, w = grid.data.shape
    stack = []
    for it in items:
        dst = np.full((h, w), np.nan, dtype="float32")
        with rasterio.Env(**GDAL_ENV), rasterio.open(signed_href(it.assets[pol].href)) as src:
            reproject(source=rasterio.band(src, 1), destination=dst, src_transform=src.transform,
                      src_crs=src.crs, src_nodata=src.nodata, dst_transform=grid.transform,
                      dst_crs=grid.crs, dst_nodata=np.nan, resampling=Resampling[resampling])
        stack.append(dst)
    if not stack:
        raise ValueError("no scene to composite")
    func = {"median": np.nanmedian, "mean": np.nanmean}[stat]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN pixels stay NaN
        out = func(np.stack(stack), axis=0).astype("float32")
    return Raster(out, grid.transform, grid.crs)


def read_cog(href: str, bbox_wgs84=None, grid: Raster | None = None,
             resampling: str = "average", nodata: float = _NODATA) -> Raster:
    """Read band 1 of a (possibly remote) GeoTIFF as linear values.

    With ``grid``, the image is resampled onto that grid (e.g. the NISAR grid, with ``average``
    in linear power, which is the correct way to aggregate backscatter); otherwise the native
    window covering ``bbox_wgs84`` is returned.
    """
    with rasterio.Env(**GDAL_ENV), rasterio.open(href) as src:
        if grid is not None:
            h, w = grid.data.shape
            with WarpedVRT(src, crs=grid.crs, transform=grid.transform, width=w, height=h,
                           resampling=Resampling[resampling], src_nodata=nodata,
                           nodata=np.nan, dtype="float32") as vrt:
                a = vrt.read(1)
            tr, crs = grid.transform, grid.crs
        else:
            if bbox_wgs84 is None:
                raise ValueError("give bbox_wgs84 or grid")
            b = bbox_to_crs(bbox_wgs84, src.crs)
            win = from_bounds(*b, transform=src.transform).round_offsets().round_lengths()
            a = src.read(1, window=win, boundless=True, fill_value=nodata).astype("float32")
            tr, crs = src.window_transform(win), src.crs
    a = a.astype("float32")
    a[(a == nodata) | (a <= 0)] = np.nan
    return Raster(a, tr, crs)


def s1_rtc_backscatter(item, bbox_wgs84=None, grid: Raster | None = None,
                       pols=("vv", "vh")) -> dict[str, Raster]:
    """Sentinel-1 gamma0 RTC in dB for each polarisation of ``item`` (keys ``"VV"``, ``"VH"``).

    Give ``grid`` to resample onto a common grid (recommended for comparisons with NISAR), or
    ``bbox_wgs84`` to read the native 10 m window.
    """
    out = {}
    for p in pols:
        if p not in item.assets:
            continue
        r = read_cog(signed_href(item.assets[p].href), bbox_wgs84=bbox_wgs84, grid=grid)
        r.data = to_db(r.data)
        out[p.upper()] = r
    if not out:
        raise KeyError(f"none of {pols} in {item.id}")
    return out
