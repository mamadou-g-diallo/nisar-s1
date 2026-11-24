# Changelog

All notable changes are listed here. Versions follow [Semantic Versioning](https://semver.org).

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
