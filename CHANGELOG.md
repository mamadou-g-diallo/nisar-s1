# Changelog

All notable changes are listed here. Versions follow [Semantic Versioning](https://semver.org).

## [0.5.0] - 2026-11-26

### Added
- `cband.common_orbit`: relative orbit present in the largest number of periods (years, seasons), so that
  multi-year Sentinel-1 comparisons keep the same viewing geometry.
- `cband.rtc_composite`: per-pixel median (or mean) of several Sentinel-1 RTC scenes on a given grid, in
  linear power. Used by the spruce budworm defoliation study (sar-tbe-defoliation), whose yearly
  composites it reproduces bit for bit.
- `grid.make_grid(..., snap=False)`: grid starting exactly at the upper-left corner of the bounds.
- Tests of these helpers (no network) and of `load_kit` with the `training/` layout.

### Changed
- `config.load_kit` looks for `training/notebooks/kit_config.yaml` first, then `notebooks/kit_config.yaml`
  (layout of the reorganised forest case-study repository).

## [0.4.2] - 2026-11-25

### Added
- `snap.s1_coherence_two_step`: Sentinel-1 coherence in two gpt runs (radar-geometry coherence cropped
  after deburst, then terrain correction), for machines with little memory. On an 8 GB laptop one
  pair took about 7 minutes instead of several hours with the single long graph. `pols="VV"` computes
  the co-polarised coherence only (the counterpart of the NISAR GUNW HH coherence).
- SNAP graphs `s1_coherence_radar` and `s1_tc`; tests of all coherence graph templates.

## [0.4.1] - 2026-11-24

### Fixed
- `workflow.read_index` keeps each file once. Indexes written by older versions could list
  the same file twice (absolute path from a first run, relative path from a second run in the same
  folder); that file then counted twice in `temporal_medians` and the `[n=...]` labels.
  Results of folders extracted in a single run are unchanged.

## [0.4.0] - 2026-11-24

First release as a standalone repository (previously developed inside the case-study repository
`nisar-s1-forest-feasibility`).

### Added
- Apache-2.0 license, NOTICE, citation file, contributing guide.
- `nisar-s1 --sites-file` and the `NISAR_S1_SITES` variable: user areas in a YAML file, on top of
  the built-in examples (`examples/sites.yaml`).
- Continuous integration with lint (ruff) and tests on Python 3.11 and 3.12.

### Changed
- Morvan example box set to the 20 x 20 km window of the case study.
- Explicit `strict=` in every `zip` and other lint fixes; no change in results.

## [0.3.0] - 2026-11

### Added
- `product`: generic reading of NISAR L2 products (content discovery, windows, memory guard,
  subsets, session files).
- `polarimetry`: dual-pol C2 from GSLC or SNAP, H/A/α, DpRVI, RVI, Stokes.
- `maps`; SNAP graph `s1_c2`.
- `nisar-s1` command line (`search`, `extract`, `summarize`, `run`, `sites`) and the `workflow`,
  `landcover`, `cband` and `raster` modules for any area.

## [0.2.0] - 2026-11

### Added
- Declared dependencies, configuration independent of the repository (`NISAR_S1_ROOT`),
  first continuous integration.
