"""Command line: ``nisar-s1 <command>`` (or ``python -m nisar_s1 <command>``).

Commands::

    nisar-s1 search    --bbox LON_MIN LAT_MIN LON_MAX LAT_MAX --start 2026-07-01 --end 2026-11-30
    nisar-s1 extract   --bbox ... --start ... --end ... --out DIR [--product GCOV GUNW S1]
    nisar-s1 summarize DIR [--target-class 10]
    nisar-s1 run       --bbox ... --start ... --end ... --out DIR      (the three steps)

``--site NAME`` can replace ``--bbox`` for the example sites of :data:`EXAMPLE_SITES` or for
your own sites, listed in a YAML file given by ``--sites-file`` or the ``NISAR_S1_SITES``
environment variable::

    my_site:
      bbox: [lon_min, lat_min, lon_max, lat_max]
      label: free text
NISAR access needs a free Earthdata account (``~/.netrc``); Sentinel-1 RTC and WorldCover are
read anonymously.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

from . import __version__

# Example boxes (lon_min, lat_min, lon_max, lat_max) covered by NISAR GCOV and GUNW products of
# the 2026 science phase: the French sites of the L-band / C-band feasibility study.
EXAMPLE_SITES = {
    "vosges": ([6.716, 47.885, 6.797, 47.935], "Vosges mountain forest, summer clear-cut (France)"),
    "haguenau": ([7.74, 48.83, 7.78, 48.86], "Haguenau temperate lowland forest (France)"),
    "morvan": ([4.07, 47.173, 4.338, 47.356], "Morvan forests, clear-cuts, 20 x 20 km window (France)"),
    "sechilienne": ([5.78, 45.03, 5.89, 45.085], "Séchilienne forested landslide (France)"),
}


def load_sites(path) -> dict[str, tuple[list[float], str]]:
    """Read user sites from a YAML file (``name: {bbox: [lon_min, lat_min, lon_max, lat_max],
    label: text}``); exits with a clear message on a malformed entry."""
    import yaml
    with open(path, encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}
    sites = {}
    for name, v in raw.items():
        bbox = (v or {}).get("bbox")
        if not (isinstance(bbox, list) and len(bbox) == 4):
            sys.exit(f"{path}: site {name!r} needs bbox: [lon_min, lat_min, lon_max, lat_max]")
        sites[str(name)] = ([float(b) for b in bbox], str(v.get("label", "")))
    return sites


def all_sites(sites_file=None) -> dict[str, tuple[list[float], str]]:
    """Example sites plus the user's (``sites_file`` or ``$NISAR_S1_SITES``); user entries win."""
    out = dict(EXAMPLE_SITES)
    f = sites_file or os.environ.get("NISAR_S1_SITES")
    if f:
        if not Path(f).is_file():
            sys.exit(f"sites file not found: {f}")
        out.update(load_sites(f))
    return out


def _bbox(args) -> list[float]:
    if args.site:
        sites = all_sites(getattr(args, "sites_file", None))
        if args.site not in sites:
            sys.exit(f"unknown site {args.site!r}; choose among {', '.join(sites)}")
        bbox = sites[args.site][0]
        x0, y0, x1, y1 = bbox
        if not (x0 < x1 and y0 < y1):
            sys.exit(f"site {args.site!r}: bbox must be lon_min, lat_min, lon_max, lat_max")
        return bbox
    if not args.bbox:
        sys.exit("give --bbox LON_MIN LAT_MIN LON_MAX LAT_MAX or --site NAME")
    x0, y0, x1, y1 = args.bbox
    if not (x0 < x1 and y0 < y1 and -180 <= x0 <= 180 and -90 <= y0 <= 90):
        sys.exit("--bbox must be LON_MIN LAT_MIN LON_MAX LAT_MAX in degrees (WGS84)")
    return [x0, y0, x1, y1]


def _area_args(p: argparse.ArgumentParser, dates: bool = True) -> None:
    p.add_argument("--bbox", nargs=4, type=float, metavar=("LON_MIN", "LAT_MIN", "LON_MAX", "LAT_MAX"))
    p.add_argument("--site", help=f"site name: {', '.join(EXAMPLE_SITES)} or one of --sites-file")
    p.add_argument("--sites-file", help="YAML file of your own sites (default: $NISAR_S1_SITES)")
    if dates:
        p.add_argument("--start", required=True, help="first date, YYYY-MM-DD")
        p.add_argument("--end", required=True, help="last date, YYYY-MM-DD")


def _extract_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--out", required=True, help="output folder (GeoTIFFs + files.csv)")
    p.add_argument("--product", nargs="+", default=["GCOV", "GUNW", "S1"],
                   choices=["GCOV", "GUNW", "S1"], help="what to extract (default: all three)")
    p.add_argument("--track", type=int, help="NISAR track (default: the most acquired)")
    p.add_argument("--direction", choices=["A", "D"], help="NISAR direction (with --track)")
    p.add_argument("--frequency", default="A", choices=["A", "B"], help="NISAR frequency")
    p.add_argument("--res", type=float, default=20.0, help="common grid posting in m (default 20)")
    p.add_argument("--max-products", type=int, help="limit per product type (quick test)")
    p.add_argument("--max-baseline", type=int, default=12,
                   help="longest GUNW temporal baseline kept, in days (default 12)")


def cmd_search(args) -> None:
    from .workflow import catalog_summary, nisar_catalog
    bbox = _bbox(args)
    print(f"Searching NISAR GCOV/GUNW over {bbox} ({args.start} -> {args.end})...", flush=True)
    df = nisar_catalog(bbox, args.start, args.end, min_cover=args.min_cover)
    if not len(df):
        print("No NISAR product fully covers this box. Try a smaller box, --min-cover 0.5, or "
              "other dates.")
        return
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(catalog_summary(df).to_string(index=False))
    print("bw_code: 2005 = 20 MHz (+5 MHz side band), 4005 = 40 MHz; pol_code: DH/DHDH = HH+HV, "
          "SH/SHSH = HH only")
    if args.csv:
        df.drop(columns=[c for c in df.columns if c == "geometry"], errors="ignore").to_csv(
            args.csv, index=False)
        print(f"Catalogue written to {args.csv}")


def _extract(args, bbox) -> None:
    from .workflow import extract_nisar, extract_s1, nisar_catalog, shared_geometry
    products = [p for p in ("GCOV", "GUNW") if p in args.product]
    if products:
        cat = nisar_catalog(bbox, args.start, args.end, products=products)
        track, direction = args.track, args.direction
        if (track is None or direction is None) and len(cat):
            track, direction = shared_geometry(cat, products)
        for product in products:
            try:
                extract_nisar(bbox, args.start, args.end, args.out, product=product,
                              track=track, direction=direction,
                              frequency=args.frequency, res=args.res, catalog=cat,
                              max_products=args.max_products,
                              max_baseline_days=args.max_baseline)
            except (ValueError, RuntimeError) as e:
                print(f"NISAR {product}: {e}")
    if "S1" in args.product:
        try:
            extract_s1(bbox, args.start, args.end, args.out, res=args.res,
                       max_products=args.max_products)
        except (ValueError, RuntimeError) as e:
            print(f"Sentinel-1: {e}")


def cmd_extract(args) -> None:
    _extract(args, _bbox(args))


def cmd_summarize(args) -> None:
    from .workflow import summarize
    stats, sep = summarize(args.out, target_class=args.target_class)
    _print_summary(stats, sep)


def _print_summary(stats: pd.DataFrame, sep: pd.DataFrame) -> None:
    if len(stats):
        piv = stats.pivot_table(index="class_name", columns="variable", values="median")
        with pd.option_context("display.width", 200, "display.max_columns", 20,
                               "display.float_format", "{:.2f}".format):
            print("\nMedian per class:\n" + piv.to_string())
    if len(sep):
        piv = sep.pivot_table(index="class_b_name", columns="variable", values="auc_sep")
        with pd.option_context("display.width", 200, "display.max_columns", 20,
                               "display.float_format", "{:.2f}".format):
            print(f"\nSeparation AUC of '{sep['class_a_name'].iloc[0]}' against each class "
                  "(0.5 = none, 1 = perfect):\n" + piv.to_string())


def cmd_run(args) -> None:
    from .workflow import summarize
    bbox = _bbox(args)
    _extract(args, bbox)
    stats, sep = summarize(args.out, target_class=args.target_class)
    _print_summary(stats, sep)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="nisar-s1", description=(
        "L band (NISAR) versus C band (Sentinel-1) on any site: catalogue, extraction on a "
        "common grid, per-land-cover statistics."))
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("search", help="list the NISAR products covering a box")
    _area_args(p)
    p.add_argument("--min-cover", type=float, default=0.95, help="fraction of the box covered")
    p.add_argument("--csv", help="write the full catalogue to this CSV file")
    p.set_defaults(func=cmd_search)

    p = sub.add_parser("extract", help="extract NISAR and Sentinel-1 layers on a common grid")
    _area_args(p)
    _extract_args(p)
    p.set_defaults(func=cmd_extract)

    p = sub.add_parser("summarize", help="statistics per WorldCover class of an output folder")
    p.add_argument("out", help="folder written by 'extract'")
    p.add_argument("--target-class", type=int, default=10,
                   help="WorldCover class compared with the others (10 tree cover, 95 mangroves)")
    p.set_defaults(func=cmd_summarize)

    p = sub.add_parser("run", help="extract then summarize")
    _area_args(p)
    _extract_args(p)
    p.add_argument("--target-class", type=int, default=10)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("sites", help="list the example sites and your own ($NISAR_S1_SITES)")
    p.add_argument("--sites-file", help="YAML file of your own sites")
    p.set_defaults(func=lambda a: print("\n".join(
        f"{k:12s} {v[1]:55s} {v[0]}" for k, v in all_sites(a.sites_file).items())))
    return ap


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
