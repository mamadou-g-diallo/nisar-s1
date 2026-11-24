# nisar-s1

[![tests](https://github.com/mamadou-g-diallo/nisar-s1/actions/workflows/tests.yml/badge.svg)](https://github.com/mamadou-g-diallo/nisar-s1/actions/workflows/tests.yml)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
![Python 3.11 | 3.12](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue.svg)

`nisar_s1` is a Python package to **read NISAR L2 products** (GCOV, GSLC, GUNW) and to **compare
L-band SAR (NISAR, later ESA ROSE-L) with C-band SAR (Sentinel-1)** on any area: catalogue search,
extraction on a common grid, coherence and backscatter statistics per land-cover class, dual-pol
polarimetry (C2 matrix, H/A/α, DpRVI, Stokes) and change separability (AUC, Cohen's d, block
bootstrap).

It was built for, and validated in, a feasibility study of clear-cut detection and forested-landslide
monitoring in France: see the case study
[nisar-s1-forest-feasibility](https://github.com/mamadou-g-diallo/nisar-s1-forest-feasibility)
(results, reports and 19 executed training notebooks).

## Installation

```bash
pip install "nisar_s1 @ git+https://github.com/mamadou-g-diallo/nisar-s1@v0.4.0"
```

For development:

```bash
git clone https://github.com/mamadou-g-diallo/nisar-s1
cd nisar-s1
pip install -e ".[dev]"
pytest
```

Python 3.11 or later. Dependencies (NumPy, h5py, rasterio, geopandas, asf_search, earthaccess, …)
are declared in `pyproject.toml`. GDAL-based packages are easiest to install from conda-forge.
ESA SNAP (`gpt`) is needed only for the Sentinel-1 SLC graphs of `nisar_s1.snap`.

## Data access and credentials

* **NISAR** (ASF DAAC): a free [NASA Earthdata](https://urs.earthdata.nasa.gov) account, in
  `~/.netrc` (`machine urs.earthdata.nasa.gov login <user> password <password>`) or in the
  `EARTHDATA_USERNAME` / `EARTHDATA_PASSWORD` variables. Products are read either from a local file
  or **directly from the HTTPS URL**: only the HDF5 chunks of the requested window are transferred.
* **Sentinel-1 RTC** and **ESA WorldCover**: public STAC catalogues (Microsoft Planetary Computer),
  no account needed.

The package never prints or stores credentials.

## Command line

```bash
nisar-s1 sites                                    # example areas
nisar-s1 search --site vosges --start 2026-06-01 --end 2026-11-20
nisar-s1 run --bbox 6.716 47.885 6.797 47.935 --start 2026-06-01 --end 2026-11-20 --out out/vosges
nisar-s1 summarize out/vosges                     # statistics of an existing output folder
```

`run` = `extract` (NISAR GCOV/GUNW and Sentinel-1 RTC on one UTM grid, GeoTIFFs plus `files.csv`)
then `summarize` (γ⁰ and coherence per WorldCover class, separability of a target class from the
others). `--max-products 2` gives a quick test.

Your own areas go in a YAML file, passed with `--sites-file` or the `NISAR_S1_SITES` variable
(see [`examples/sites.yaml`](examples/sites.yaml)):

```yaml
my_forest:
  bbox: [lon_min, lat_min, lon_max, lat_max]   # WGS84
  label: "Short description"
```

## Python API

| Module | Purpose |
|---|---|
| `product` | Discover a NISAR L2 product from the file itself (`product_profile`: type, frequencies, polarisations, grid), metadata without pixel reads, windowed reads with a memory guard (`read_layer`, `block_iter`), AOI subsets (`write_subset`), session files |
| `polarimetry` | Dual-pol C2 matrix from GSLC (`c2_from_slc`) or SNAP elements, physical validation, eigenvalues, H/A/ᾱ (Cloude 2007), DpRVI (Mandal et al. 2020), RVI, Stokes parameters, RGB composites |
| `io_nisar` | Open local or remote NISAR products, GCOV backscatter and GUNW coherence on a window |
| `search` | ASF search for NISAR and Sentinel-1, name parsing, revisit summary, pair chains |
| `cband` | Sentinel-1 RTC backscatter from the Planetary Computer |
| `workflow` | Catalogue, geometry choice, extraction on a reference grid, temporal medians, summary (what the CLI runs) |
| `grid`, `raster` | Common grid, reprojection, polygon masks, GeoTIFF I/O, dB conversion |
| `landcover` | WorldCover classes, erosion of class borders, per-class statistics and separability |
| `stats` | AUC, Cohen's d, Jeffries–Matusita, block bootstrap of the AUC (spatial autocorrelation) |
| `plots`, `maps` | Figures (time series, separability bars, change maps, matrices), static and folium maps |
| `snap` | Sentinel-1 SLC graphs for ESA SNAP `gpt`: γ⁰, coherence, interferogram, C2 matrix |

Example: dual-pol H/A/α from two GSLC polarisations

```python
import h5py
from nisar_s1 import product, polarimetry as pol

with h5py.File("NISAR_L2_PR_GSLC_....h5") as h5:
    prof = product.product_profile(h5)
    grid = prof.grids["A"]
    win = product.window_from_bounds(grid, (x0, y0, x1, y1))   # in the product CRS
    hh = product.read_layer(h5, grid, "HH", win)
    hv = product.read_layer(h5, grid, "HV", win)

c2 = pol.c2_from_slc(hh, hv, looks=(8, 4))
print(pol.validate_c2(c2))          # fractions of non-physical pixels, should be ~0
haa = pol.h_alpha_dual(c2)          # dict of arrays: entropy, anisotropy, mean alpha
```

The case-study helpers (`config`, `download`, `optical`, `sites`, `ign`, `detect`, `guide_cr`)
implement the French clear-cut workflow (Sentinel-2 reference, IGN data, the Irstea clear-cut
detection guide of Ose & Deshayes 2015). They are kept public but are less general.

## Tests

```bash
pytest
ruff check src tests
```

The tests use synthetic products (small HDF5 files with the NISAR layout) and need no network.
Continuous integration runs them on Python 3.11 and 3.12.

## Citation

If you use this package, please cite it (see [`CITATION.cff`](CITATION.cff)):

> Diallo, M. G. (2026). *nisar-s1: NISAR L2 reading and L-band / C-band SAR comparison tools*
> (version 0.4.0). https://github.com/mamadou-g-diallo/nisar-s1

## Acknowledgements

* NISAR L2 products are generated by NASA JPL with **ISCE3**
  ([isce-framework/isce3](https://github.com/isce-framework/isce3)) and distributed by the ASF DAAC.
  NISAR is a joint mission of NASA and ISRO.
* The package is designed so that the same code applies to the L2 products of ESA's future L-band
  mission **ROSE-L**.
* The structure of the training notebooks of the case study follows the NISAR GCOV training kit of
  D. K. Pandey ([dharmendraisro/nisar-polarimetric-training](https://github.com/dharmendraisro/nisar-polarimetric-training),
  Apache-2.0); no code is copied from it.
* Sentinel-1 and Sentinel-2 data: Copernicus programme (ESA / European Union). ESA WorldCover.

## License

Apache License 2.0, see [LICENSE](LICENSE) and [NOTICE](NOTICE).
