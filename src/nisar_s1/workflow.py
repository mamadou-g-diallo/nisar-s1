"""High-level workflow for a new site: catalogue, extraction on a common grid, class summary.

Three steps, each usable from Python or from the ``nisar-s1`` command line (:mod:`nisar_s1.cli`):

1. :func:`nisar_catalog` - which NISAR products (GCOV backscatter, GUNW interferometry) cover the
   box, in which geometry (track / direction) and acquisition mode;
2. :func:`extract_nisar` and :func:`extract_s1` - read only the box from each product (remote
   HTTP range requests, nothing large is downloaded) and write GeoTIFFs on one common grid, plus
   an index ``files.csv``;
3. :func:`summarize` - temporal median of each variable, WorldCover classes on the same grid,
   per-class statistics and the separation (AUC) of one class against the others, L band versus
   C band.

The grid of the first NISAR product read defines the common grid (its UTM zone and posting);
every later product, and Sentinel-1, is resampled onto it (average of linear power for
backscatter), so per-pixel comparisons use the same ground samples.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd

from .config import bbox_to_wkt
from .grid import make_grid, to_grid
from .raster import Raster, bbox_to_crs, load_raster, save_raster, to_db

INDEX = "files.csv"


def _say(msg: str) -> None:
    print(msg, flush=True)


# ---------------------------------------------------------------------------------- catalogue
def nisar_catalog(bbox_wgs84, start: str, end: str, products=("GCOV", "GUNW"),
                  min_cover: float = 0.95) -> pd.DataFrame:
    """NISAR products fully (``min_cover``) covering the box, with geometry and mode columns."""
    from .search import search_nisar
    df = search_nisar(bbox_to_wkt(bbox_wgs84), start, end, product_types=tuple(products))
    if not len(df):
        return df
    df = df[df["aoi_cover"] >= min_cover].copy()
    df["geometry_id"] = df["track"].astype(str) + df["direction"]
    return df.sort_values("startTime").reset_index(drop=True)


def catalog_summary(df: pd.DataFrame) -> pd.DataFrame:
    """One row per product / geometry / mode: number of products and date range."""
    if not len(df):
        return pd.DataFrame()
    g = df.groupby(["product", "track", "direction", "frame", "bw_code", "pol_code"])
    out = g.agg(n=("sceneName", "size"), first=("ref_date", "min"), last=("ref_date", "max"))
    if "baseline_days" in df:
        out["baselines_d"] = g["baseline_days"].agg(
            lambda s: ",".join(str(int(v)) for v in sorted(s.dropna().unique())))
    return out.reset_index()


def best_geometry(df: pd.DataFrame, product: str = "GCOV") -> tuple[int, str]:
    """(track, direction) with the most acquisitions of ``product`` (dual-pol preferred)."""
    d = df[df["product"] == product]
    if not len(d):
        raise ValueError(f"no {product} product in the catalogue")
    score = d.groupby(["track", "direction"]).apply(
        lambda g: (g["pol_code"].str.startswith("D").sum(), len(g)), include_groups=False)
    track, direction = max(score.items(), key=lambda kv: kv[1])[0]
    return int(track), str(direction)


def shared_geometry(df: pd.DataFrame, products=("GCOV", "GUNW")) -> tuple[int, str]:
    """One (track, direction) for all ``products``, so backscatter and coherence share the
    viewing geometry: the geometry with the most acquisitions among those that have every
    product (dual-pol GCOV preferred), else the best GCOV (or only) geometry."""
    products = [p for p in products if p in set(df["product"])]
    if not products:
        raise ValueError("empty catalogue")
    have = set.intersection(*(set(zip(df[df["product"] == p]["track"],
                                      df[df["product"] == p]["direction"], strict=True)) for p in products))
    main = "GCOV" if "GCOV" in products else products[0]
    d = df[df["product"] == main]
    if have:
        d = d[[(t, r) in have for t, r in zip(d["track"], d["direction"], strict=True)]]
    return best_geometry(d, main)


# --------------------------------------------------------------------------------- extraction
def read_index(out_dir) -> pd.DataFrame:
    """The index of an output folder, with absolute paths.

    Paths are stored as file names relative to the folder, so the folder can be moved or read
    from another working directory. Older indexes (absolute paths, or paths relative to the
    directory the extraction was run from) are resolved the same way. A file listed twice (older
    indexes mixing absolute and relative paths after a second run) is kept once, with its last
    entry, so that no date counts twice in the temporal medians.
    """
    out_dir = Path(out_dir).resolve()
    idx = pd.read_csv(out_dir / INDEX)
    idx["path"] = [str(p) if Path(p).is_absolute() and Path(p).exists() else str(out_dir / Path(p).name)
                   for p in idx["path"]]
    name = idx["path"].map(lambda p: Path(p).name)
    return idx[~name.duplicated(keep="last")].reset_index(drop=True)


def _append_index(out_dir: Path, rows: list[dict]) -> pd.DataFrame:
    path = out_dir / INDEX
    new = pd.DataFrame(rows)
    new["path"] = [Path(p).name for p in new["path"]]
    if path.exists():
        old = pd.read_csv(path)
        old["path"] = [Path(p).name for p in old["path"]]
        new = pd.concat([old[~old["path"].isin(new["path"])], new], ignore_index=True)
    new.to_csv(path, index=False)
    return read_index(out_dir)


def _regrid(r: Raster, grid: Raster, db: bool) -> Raster:
    """Resample onto ``grid`` (pixel average). Backscatter is averaged in linear power, not in dB."""
    same = (r.data.shape == grid.data.shape and r.transform == grid.transform
            and str(r.crs) == str(grid.crs))
    if same:
        return r
    if not db:
        return to_grid(r, grid)
    lin = Raster((10 ** (r.data.astype("float64") / 10)).astype("float32"), r.transform, r.crs)
    out = to_grid(lin, grid)
    out.data = to_db(out.data)
    return out


def extract_nisar(bbox_wgs84, start: str, end: str, out_dir, product: str = "GCOV",
                  track: int | None = None, direction: str | None = None,
                  frequency: str = "A", kind: str = "gamma0", coherence_posting: int = 20,
                  max_baseline_days: int | None = 12, max_products: int | None = None,
                  res: float = 20.0, catalog: pd.DataFrame | None = None,
                  grid: Raster | None = None,
                  open_fn: Callable | None = None, say: Callable = _say) -> pd.DataFrame:
    """Read the box from every NISAR ``product`` (``GCOV`` or ``GUNW``) of one geometry and write
    GeoTIFFs in ``out_dir`` (backscatter in dB per polarisation, or coherence 0-1).

    The geometry defaults to :func:`best_geometry`. GUNW pairs longer than
    ``max_baseline_days`` are skipped (coherence is only comparable at equal baselines).
    Every layer is resampled onto the common grid (``grid``, else the grid already used in
    ``out_dir``, else ``res`` m in the product projection; 20 m by default, the GUNW coherence
    posting, which also averages 2x2 looks of the 10 m GCOV of the 40 MHz mode).
    Existing files are kept (the extraction can be resumed). Returns the index of written files.
    """
    from .io_nisar import gcov_backscatter, gunw_coherence, open_product
    open_fn = open_fn or open_product
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    grid = grid if grid is not None else reference_grid(out_dir)
    df = catalog if catalog is not None else nisar_catalog(bbox_wgs84, start, end, (product,))
    df = df[df["product"] == product]
    if not len(df):
        raise ValueError(f"no NISAR {product} product covers {list(bbox_wgs84)} "
                         f"between {start} and {end}")
    if track is None or direction is None:
        track, direction = best_geometry(df, product)
    df = df[(df["track"] == track) & (df["direction"] == direction)]
    if product == "GUNW" and max_baseline_days:
        df = df[df["baseline_days"] <= max_baseline_days]
    df = df.drop_duplicates("sceneName").sort_values("startTime")
    if max_products:
        df = df.head(max_products)
    geom = f"T{track:03d}{direction}"
    say(f"NISAR {product} {geom}: {len(df)} product(s) to read")
    rows = []
    for k, (_, p) in enumerate(df.iterrows(), 1):
        d1 = pd.Timestamp(p["ref_date"]).strftime("%Y%m%d")
        d2 = pd.Timestamp(p["sec_date"]).strftime("%Y%m%d") if product == "GUNW" else ""
        stem = f"nisar_{product}_{geom}_{d1}" + (f"_{d2}" if d2 else "")
        say(f"  [{k}/{len(df)}] {stem}")
        if product == "GCOV":
            targets = {pol: out_dir / f"{stem}_{pol}_{kind}_dB.tif" for pol in ("HH", "HV", "VV", "VH")}
        else:
            targets = {f"coh{coherence_posting}": out_dir / f"{stem}_coh{coherence_posting}.tif"}
        layers = {v: load_raster(t) for v, t in targets.items() if t.exists()}
        if layers:
            say("      already extracted")
        else:
            try:
                with open_fn(p["url"]) as h5:
                    if product == "GCOV":
                        layers = gcov_backscatter(h5, kind=kind, bbox_wgs84=bbox_wgs84,
                                                  frequency=frequency)
                    else:
                        layers = {f"coh{coherence_posting}": gunw_coherence(
                            h5, bbox_wgs84=bbox_wgs84, posting_m=coherence_posting,
                            frequency=frequency)}
            except (OSError, KeyError, ValueError) as e:
                say(f"      skipped: {e}")
                continue
        for var, r in layers.items():
            if grid is None:
                grid = make_grid(bbox_to_crs(bbox_wgs84, r.crs), res, str(r.crs))
            path = targets[var]
            if not path.exists():
                save_raster(_regrid(r, grid, db=product == "GCOV"), path,
                            description=f"NISAR {product} {var}")
            rows.append({"sensor": "NISAR", "band": "L", "product": product, "geometry": geom,
                         "date": pd.Timestamp(p["ref_date"]).date(),
                         "date2": pd.Timestamp(p["sec_date"]).date() if d2 else None,
                         "variable": f"L {var}" + (" (dB)" if product == "GCOV" else ""),
                         "path": str(path)})
    if not rows:
        raise RuntimeError("no product could be read")
    return _append_index(out_dir, rows)


def reference_grid(out_dir) -> Raster | None:
    """Grid of the first file of the index (``None`` if nothing was extracted yet)."""
    path = Path(out_dir) / INDEX
    if not path.exists():
        return None
    first = read_index(out_dir)["path"].iloc[0]
    r = load_raster(first)
    return Raster(np.full(r.data.shape, np.nan, "float32"), r.transform, r.crs)


def extract_s1(bbox_wgs84, start: str, end: str, out_dir, orbit_state: str | None = None,
               relative_orbit: int | None = None, grid: Raster | None = None,
               res: float = 20.0, max_products: int | None = None,
               say: Callable = _say) -> pd.DataFrame:
    """Sentinel-1 gamma0 RTC (VV, VH, dB) for every scene of one relative orbit, resampled onto
    the common grid (the NISAR grid already in ``out_dir`` if any, else ``res`` m in UTM)."""
    from .cband import s1_rtc_backscatter, search_s1_rtc
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = search_s1_rtc(bbox_wgs84, start, end, orbit_state=orbit_state,
                       relative_orbit=relative_orbit)
    if not len(df):
        raise ValueError(f"no Sentinel-1 RTC scene covers {list(bbox_wgs84)} "
                         f"between {start} and {end}")
    if relative_orbit is None:
        relative_orbit = int(df["relative_orbit"].value_counts().idxmax())
        df = df[df["relative_orbit"] == relative_orbit]
    df = df.assign(day=df["datetime"].dt.date).drop_duplicates("day")
    if max_products:
        df = df.head(max_products)
    grid = grid if grid is not None else reference_grid(out_dir)
    if grid is None:
        from .cband import read_cog, signed_href
        item = df.iloc[0]["item"]
        r = read_cog(signed_href(item.assets["vv"].href), bbox_wgs84=bbox_wgs84)
        grid = make_grid(bbox_to_crs(bbox_wgs84, r.crs), res, str(r.crs))
    geom = f"R{relative_orbit:03d}{str(df.iloc[0]['orbit_state'])[:1].upper()}"
    say(f"Sentinel-1 RTC {geom}: {len(df)} scene(s) to read")
    rows = []
    for k, (_, s) in enumerate(df.iterrows(), 1):
        d = s["datetime"].strftime("%Y%m%d")
        say(f"  [{k}/{len(df)}] s1_rtc_{geom}_{d}")
        targets = {p: out_dir / f"s1_RTC_{geom}_{d}_{p}_gamma0_dB.tif" for p in ("VV", "VH")}
        try:
            missing = [p.lower() for p, t in targets.items() if not t.exists()]
            layers = s1_rtc_backscatter(s["item"], grid=grid, pols=missing) if missing else {}
        except Exception as e:  # noqa: BLE001 - network or token errors: skip the scene
            import re
            say(f"      skipped: {re.sub(r'[?][^ ]*', '?<token>', str(e))[:300]}")
            continue
        for p, path in targets.items():
            if p in layers:
                save_raster(layers[p], path, description=f"Sentinel-1 RTC {p}")
            if path.exists():
                rows.append({"sensor": "Sentinel-1", "band": "C", "product": "RTC",
                             "geometry": geom, "date": s["datetime"].date(), "date2": None,
                             "variable": f"C {p} (dB)", "path": str(path)})
    if not rows:
        raise RuntimeError("no scene could be read")
    return _append_index(out_dir, rows)


# ------------------------------------------------------------------------------------ summary
def temporal_medians(out_dir) -> dict[str, Raster]:
    """Per-pixel temporal median of every variable of the index (one Raster per variable)."""
    idx = read_index(out_dir)
    out = {}
    for var, g in idx.groupby("variable", sort=False):
        rs = [load_raster(p) for p in g["path"]]
        shapes = {r.data.shape for r in rs}
        if len(shapes) > 1:
            raise ValueError(f"{var}: files are not on one grid {shapes}")
        with np.errstate(all="ignore"):
            med = np.nanmedian(np.stack([r.data for r in rs]), axis=0).astype("float32")
        out[f"{var} [n={len(rs)}]"] = Raster(med, rs[0].transform, rs[0].crs)
    return out


def summarize(out_dir, target_class: int = 10, erode_pixels: int = 1,
              classes: np.ndarray | None = None, figure: bool = True,
              say: Callable = _say) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Class statistics and separability of ``target_class`` (default WorldCover tree cover)
    against every other class, for each variable; writes ``class_stats.csv``,
    ``separability.csv`` and ``class_profiles.png`` in ``out_dir``."""
    from .landcover import WORLDCOVER_CLASSES, class_separability, class_stats, erode_classes, worldcover_like
    out_dir = Path(out_dir)
    layers = temporal_medians(out_dir)
    ref = next(iter(layers.values()))
    if classes is None:
        say("WorldCover classes on the common grid...")
        classes = worldcover_like(ref)
    if erode_pixels:
        classes = erode_classes(classes, erode_pixels)
    stats = class_stats(layers, classes)
    stats.to_csv(out_dir / "class_stats.csv", index=False)
    seps = []
    present = [c for c in np.unique(classes) if c and c != target_class
               and (classes == c).sum() >= 50]
    for c in present:
        seps.append(class_separability(layers, classes, target_class, int(c)))
    sep = pd.concat(seps, ignore_index=True) if seps else pd.DataFrame()
    if len(sep):
        sep["class_a_name"] = sep["class_a"].map(WORLDCOVER_CLASSES)
        sep["class_b_name"] = sep["class_b"].map(WORLDCOVER_CLASSES)
    sep.to_csv(out_dir / "separability.csv", index=False)
    if figure and len(stats):
        from .plots import class_profiles
        class_profiles(stats, out_dir / "class_profiles.png")
    say(f"Written: {out_dir / 'class_stats.csv'}, {out_dir / 'separability.csv'}"
        + (f", {out_dir / 'class_profiles.png'}" if figure else ""))
    return stats, sep
