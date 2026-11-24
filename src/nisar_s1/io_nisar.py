"""Readers for NISAR L2 geocoded products (GCOV backscatter, GUNW interferometry).

Products are HDF5 files on a projected grid (UTM) with 512x512 chunks, so a small window can be
read from a local file or directly from the Earthdata HTTPS URL (only the needed chunks are
transferred). The same API is meant to be reused for ROSE-L L2 products, which are expected to
follow a similar layout (geocoded grids + coordinate vectors + projection EPSG).

What varies between acquisitions, and is handled here:

- **Frequency**: ``A`` is the main band (20 MHz mode -> 20 m GCOV posting; 40 MHz -> 10 m),
  ``B`` the 5 MHz side band (80 m GCOV) when acquired (``listOfFrequencies``).
- **Polarisations**: dual-pol ``HH/HV`` or single-pol ``HH`` (and other modes elsewhere), read
  from ``listOfPolarizations`` / ``listOfCovarianceTerms`` instead of being assumed.
- **GUNW layers**: ``wrappedInterferogram/<pol>/coherenceMagnitude`` (20 m) and
  ``unwrappedInterferogram/<pol>/{coherenceMagnitude,unwrappedPhase,connectedComponents,
  ionospherePhaseScreen}`` (80 m).

Windows can be given in the product projection (``bounds``) or as a WGS84 box (``bbox_wgs84``).
"""

from __future__ import annotations

import contextlib
import zlib
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
from affine import Affine

from .raster import Raster, bbox_to_crs, to_db
from .search import parse_nisar_name

GCOV_GRID = "science/LSAR/GCOV/grids/frequencyA"
GUNW_GRID = "science/LSAR/GUNW/grids/frequencyA"

# short name -> (group relative to the product grid, dataset); "{pol}" is filled in by the reader
LAYERS = {
    "GCOV": {
        "HHHH": ("", "HHHH"), "HVHV": ("", "HVHV"), "VVVV": ("", "VVVV"), "VHVH": ("", "VHVH"),
        "gamma_to_sigma": ("", "rtcGammaToSigmaFactor"),
        "looks": ("", "numberOfLooks"), "mask": ("", "mask"),
    },
    "GUNW": {
        "coh20": ("wrappedInterferogram", "{pol}/coherenceMagnitude"),
        "coh80": ("unwrappedInterferogram", "{pol}/coherenceMagnitude"),
        "unw": ("unwrappedInterferogram", "{pol}/unwrappedPhase"),
        "conncomp": ("unwrappedInterferogram", "{pol}/connectedComponents"),
        "iono": ("unwrappedInterferogram", "{pol}/ionospherePhaseScreen"),
        "mask20": ("wrappedInterferogram", "mask"),
    },
}


@dataclass
class NisarInfo:
    """Acquisition metadata parsed from the granule name."""
    name: str
    product: str
    track: int
    direction: str
    frame: int
    ref_date: datetime
    sec_date: datetime | None = None

    @property
    def baseline_days(self) -> int | None:
        return None if self.sec_date is None else (self.sec_date - self.ref_date).days


def info_from_name(name: str) -> NisarInfo:
    p = parse_nisar_name(Path(name).name)
    return NisarInfo(Path(name).stem, p["product"], p["track"], p["direction"], p["frame"],
                     p["ref_date"], p.get("sec_date"))


@contextlib.contextmanager
def open_product(src: str | Path, fs=None, block_size: int = 2**20):
    """Open a NISAR HDF5 product from a local path or an HTTPS URL (``fs``: an authenticated
    fsspec filesystem, e.g. ``earthaccess.get_fsspec_https_session()``; by default Earthdata
    credentials are read from ``~/.netrc`` or the ``EARTHDATA_USERNAME/PASSWORD`` variables)."""
    if str(src).startswith("http"):
        if fs is None:
            fs = earthdata_fs()
        fobj = fs.open(str(src), block_size=block_size, cache_type="blockcache")
        h5 = h5py.File(fobj, "r")
        _REMOTE[h5.id.id] = (fs, str(src))
        try:
            yield h5
        finally:
            _REMOTE.pop(h5.id.id, None)
            h5.close()
            fobj.close()
    else:
        with h5py.File(src, "r") as h5:
            yield h5


_FS = None
# open remote products: HDF5 file id -> (fsspec filesystem, URL), used by the parallel reader
_REMOTE: dict[int, tuple] = {}
_FILTERS_OK = {"shuffle", "gzip", "fletcher32"}


def _read_chunks_parallel(dset: h5py.Dataset, i0: int, i1: int, j0: int, j1: int):
    """Read ``dset[i0:i1, j0:j1]`` of a remote product by fetching only the compressed chunks of
    the window, all at once (concurrent HTTP range requests), and decoding them here.

    Reading through h5py on a remote file fetches the chunks one after the other through a
    block cache, which is 50 to 100 times slower for a window of tens of chunks. Returns
    ``None`` when the fast path does not apply (local file, unsupported filter), so the caller
    falls back to h5py; the result is identical to ``dset[i0:i1, j0:j1]`` (tested).
    """
    remote = _REMOTE.get(dset.file.id.id)
    filters = set(getattr(dset, "_filters", {}) or {})
    if remote is None or dset.chunks is None or dset.ndim != 2 or not filters <= _FILTERS_OK:
        return None
    fs, url = remote
    ch = dset.chunks
    coords = [(ci, cj) for ci in range(i0 // ch[0] * ch[0], i1, ch[0])
              for cj in range(j0 // ch[1] * ch[1], j1, ch[1])]
    infos = [dset.id.get_chunk_info_by_coord(c) for c in coords]
    if any(i.byte_offset is not None and i.filter_mask for i in infos):
        return None  # a filter was skipped on some chunk: let HDF5 decode
    stored = [(c, i) for c, i in zip(coords, infos, strict=True) if i.byte_offset is not None]
    blobs = fs.cat_ranges([url] * len(stored), [i.byte_offset for _, i in stored],
                          [i.byte_offset + i.size for _, i in stored])
    fill = dset.fillvalue if dset.fillvalue is not None else 0
    out = np.full((i1 - i0, j1 - j0), fill, dset.dtype)
    n = dset.dtype.itemsize
    for ((ci, cj), _), raw in zip(stored, blobs, strict=True):
        if "fletcher32" in filters:
            raw = raw[:-4]
        b = zlib.decompress(raw) if "gzip" in filters else raw
        if "shuffle" in filters and n > 1:
            b = np.frombuffer(b, np.uint8).reshape(n, -1).T.tobytes()
        a = np.frombuffer(b, dset.dtype).reshape(ch)
        r0, c0 = max(ci, i0), max(cj, j0)
        r1, c1 = min(ci + ch[0], i1), min(cj + ch[1], j1)
        out[r0 - i0:r1 - i0, c0 - j0:c1 - j0] = a[r0 - ci:r1 - ci, c0 - cj:c1 - cj]
    return out


def earthdata_fs():
    """Authenticated fsspec HTTPS filesystem for Earthdata (created once per session)."""
    global _FS
    if _FS is None:
        import earthaccess
        auth = None
        for strategy in ("netrc", "environment"):
            try:
                auth = earthaccess.login(strategy=strategy)
            except Exception:  # noqa: BLE001 - try the next strategy
                auth = None
            if auth is not None and getattr(auth, "authenticated", False):
                break
        else:
            raise RuntimeError(
                "Earthdata login failed: put 'machine urs.earthdata.nasa.gov login <user> "
                "password <pwd>' in ~/.netrc, or set EARTHDATA_USERNAME and EARTHDATA_PASSWORD. "
                "A free account is created at https://urs.earthdata.nasa.gov.")
        _FS = earthaccess.get_fsspec_https_session()
    return _FS


def _root(product: str, frequency: str = "A") -> str:
    return f"science/LSAR/{product}/grids/frequency{frequency}"


def _grid(h5: h5py.File, product: str, group: str, frequency: str = "A") -> h5py.Group:
    root = _root(product, frequency)
    if root not in h5:
        raise KeyError(f"{product} frequency {frequency} not in this product "
                       f"(frequencies: {frequencies(h5, product)})")
    return h5[f"{root}/{group}" if group else root]


def _strings(a) -> list[str]:
    return [v.decode() if isinstance(v, bytes) else str(v) for v in np.atleast_1d(a)]


def frequencies(h5: h5py.File, product: str | None = None) -> list[str]:
    """Frequencies present in the product (``["A"]`` or ``["A", "B"]``)."""
    key = "science/LSAR/identification/listOfFrequencies"
    if key in h5:
        return _strings(h5[key][()])
    product = product or product_type(h5)
    grids = h5.get(f"science/LSAR/{product}/grids", {})
    return sorted(k.replace("frequency", "") for k in grids)


def product_type(h5: h5py.File) -> str:
    """``GCOV`` or ``GUNW`` (from the HDF5 layout)."""
    for p in ("GCOV", "GUNW"):
        if f"science/LSAR/{p}" in h5:
            return p
    raise ValueError("not a NISAR GCOV or GUNW product")


def polarizations(h5: h5py.File, product: str, frequency: str = "A") -> list[str]:
    """Polarisations available for ``frequency`` (e.g. ``["HH", "HV"]`` or ``["HH"]``)."""
    g = _grid(h5, product, "", frequency)
    if "listOfPolarizations" in g:
        return _strings(g["listOfPolarizations"][()])
    if product == "GCOV":
        return [k[:2] for k in ("HHHH", "HVHV", "VVVV", "VHVH") if k in g]
    sub = g.get("wrappedInterferogram", g.get("unwrappedInterferogram", {}))
    return [k for k in ("HH", "HV", "VV", "VH") if k in sub]


def posting(h5: h5py.File, product: str, layer: str = "", frequency: str = "A") -> float:
    """Pixel spacing (m) of a GCOV frequency, or of a GUNW layer (``coh20``/``coh80``...)."""
    group = LAYERS[product][layer][0] if layer else ""
    g = _grid(h5, product, group, frequency)
    x = g["xCoordinates"]
    return float(abs(x[1] - x[0]))


def epsg(h5: h5py.File, product: str | None = None, frequency: str = "A") -> int:
    """EPSG code of the product grid."""
    product = product or product_type(h5)
    g = _grid(h5, product, "" if product == "GCOV" else "wrappedInterferogram", frequency)
    return int(g["projection"][()])


def read_window(h5: h5py.File, product: str, layer: str, bounds=None, *, bbox_wgs84=None,
                frequency: str = "A", pol: str = "HH") -> Raster:
    """Read ``layer`` (a key of :data:`LAYERS`) inside a window.

    The window is ``bounds`` = (xmin, ymin, xmax, ymax) in the product projection, or
    ``bbox_wgs84`` = [lon_min, lat_min, lon_max, lat_max]. Returns a :class:`Raster` (NaN where
    float layers have no data; NISAR uses 0 as the no-data value of float layers).
    """
    group, dset = LAYERS[product][layer]
    dset = dset.format(pol=pol)
    g = _grid(h5, product, group, frequency)
    if dset not in g:
        raise KeyError(f"{product} layer {dset!r} not found for frequency {frequency} "
                       f"(polarisations: {polarizations(h5, product, frequency)})")
    x, y = g["xCoordinates"][:], g["yCoordinates"][:]
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])  # dy < 0
    crs = f"EPSG:{int(g['projection'][()])}"
    if bounds is None:
        if bbox_wgs84 is None:
            raise ValueError("give bounds (product projection) or bbox_wgs84")
        bounds = bbox_to_crs(bbox_wgs84, crs)
    xmin, ymin, xmax, ymax = bounds
    j0 = max(int(np.floor((xmin - x[0]) / dx + 0.5)), 0)
    j1 = min(int(np.ceil((xmax - x[0]) / dx + 0.5)), len(x))
    i0 = max(int(np.floor((ymax - y[0]) / dy + 0.5)), 0)
    i1 = min(int(np.ceil((ymin - y[0]) / dy + 0.5)), len(y))
    if j1 <= j0 or i1 <= i0:
        raise ValueError(f"window {tuple(round(b) for b in bounds)} outside the product grid")
    data = _read_chunks_parallel(g[dset], i0, i1, j0, j1)
    if data is None:
        data = g[dset][i0:i1, j0:j1]
    if np.issubdtype(data.dtype, np.floating):
        data = data.astype("float32")
        data[data == 0] = np.nan
    # coordinates are pixel centres
    transform = Affine(dx, 0, x[j0] - dx / 2, 0, dy, y[i0] - dy / 2)
    return Raster(data, transform, crs)


def gcov_backscatter(h5: h5py.File, bounds=None, kind: str = "gamma0", *, bbox_wgs84=None,
                     frequency: str = "A", pols=None) -> dict[str, Raster]:
    """Backscatter in dB for every available co-/cross-polarisation (e.g. ``{"HH", "HV"}``).

    ``kind`` = ``gamma0`` (radiometrically terrain corrected, the GCOV layer) or ``sigma0``
    (gamma0 x ``rtcGammaToSigmaFactor``). ``pols`` restricts the output (default: all
    diagonal covariance terms present, so single-pol products return ``{"HH"}`` only).
    """
    if kind not in ("gamma0", "sigma0"):
        raise ValueError("kind must be 'gamma0' or 'sigma0'")
    available = polarizations(h5, "GCOV", frequency)
    wanted = available if pols is None else [p for p in pols if p in available]
    if not wanted:
        raise KeyError(f"none of {pols} in this product (available: {available})")
    kw = dict(bounds=bounds, bbox_wgs84=bbox_wgs84, frequency=frequency)
    fac = read_window(h5, "GCOV", "gamma_to_sigma", **kw).data if kind == "sigma0" else 1.0
    out = {}
    for pol in wanted:
        r = read_window(h5, "GCOV", pol * 2, **kw)
        r.data = to_db(r.data * fac)
        out[pol] = r
    return out


def gunw_coherence(h5: h5py.File, bounds=None, *, bbox_wgs84=None, posting_m: int = 20,
                   frequency: str = "A", pol: str = "HH") -> Raster:
    """Interferometric coherence magnitude (0-1) of a GUNW product, at 20 m (wrapped
    interferogram layer) or 80 m (unwrapped layer, more looks hence lower bias)."""
    layer = {20: "coh20", 80: "coh80"}.get(int(posting_m))
    if layer is None:
        raise ValueError("posting_m must be 20 or 80")
    return read_window(h5, "GUNW", layer, bounds, bbox_wgs84=bbox_wgs84, frequency=frequency,
                       pol=pol)
