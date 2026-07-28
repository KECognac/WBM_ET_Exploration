# WBM_ET_Exploration

Validating the National Park Service (NPS) Water Balance Model (WBM) against
OpenET remote-sensing evapotranspiration and AmeriFlux/USGS eddy-covariance
flux towers, across the full nationwide list of flux tower sites with data
overlapping the OpenET period (2016–2023) — roughly 80 sites spanning most
major CONUS ecosystems. A separate side analysis (notebook 08) restricts the
same cached results to the smaller subset of sites within 50 km of an NPS
park unit, which was this project's original scope.

Originally developed as a class project (CIVE 523, Remote Sensing Hydrology)
around six NPS sites; reorganized here into a standalone, reproducible repo
and later expanded to the nationwide flux tower list so the same validation
can inform how the WBM is used across the NPS system generally, not just at
those six parks.

### 📊 [View the results report](https://nps-water-supply-vulnerability.github.io/WBM_ET_Exploration/reports/pet_comparison_report.html)

A concise, self-contained HTML summary comparing the Oudin and
Penman-Monteith WBM AET runs against OpenET and flux towers (live page via
GitHub Pages — see [Setup](#setup) below if it's not enabled yet). Source
file: [`reports/pet_comparison_report.html`](reports/pet_comparison_report.html).

## What this does

1. Selects the full nationwide list of AmeriFlux/USGS flux towers with data
   overlapping the OpenET period (2016–2023) — not restricted to sites near a
   park — and saves it to `Data/geo_data/flux_towers.csv`. Separately
   identifies the subset of those sites within 50 km of an NPS park unit
   (`Data/geo_data/flux_towers_near_parks.csv`), used only by notebook 08's
   side analysis below.
2. Pulls observed flux-tower ET, point-sampled and grid-cell-averaged OpenET
   (via Google Earth Engine), and gridMET daily climate for every site in the
   nationwide list.
3. Runs the NPS WBM (Oudin PET, degree-day snowmelt, bucket-style soil
   moisture) at a daily timestep, driven by gridMET.
4. Compares WBM actual ET (AET) to OpenET and flux-tower ET using bias,
   RMSE, and NSE, then attempts bias correction (direct AET correction, and
   a per-site PET-multiplier calibration).
5. Produces before/after comparison figures and exports a final metrics
   table.
6. (`notebook 06`) Runs a second, uncalibrated WBM configuration driven by
   FAO-56 Penman-Monteith instead of Oudin, and compares the two PET/AET
   estimates head-to-head. Penman-Monteith returns ETo (reference grass ET),
   not PET, and needs real wind speed — gridMET's wind (`vs` band, 10 m) is
   downloaded in notebook 02 alongside precip/temp and corrected there to
   the 2 m height FAO-56 assumes (`wind_ms` column). Notebook 06 applies a
   crop coefficient (`Kc = 1.0` baseline, easy to change) to convert
   ETo → PET before running it through the same WBM soil/snow/AET machinery
   as the Oudin run.
7. (`notebook 07`) Benchmarks *three* WBM AET variants — default Oudin,
   notebook 04's PET-multiplier-**calibrated** Oudin, and (uncalibrated)
   Penman-Monteith — against OpenET (the same 4 km gridcell ensemble
   notebook 04 validates against) and against flux tower ET, producing
   per-site bias/RMSE/R² tables for all six comparisons. At nationwide scale
   (~80 sites) the figures are condensed rather than faceted per site: one
   scatter plot per AET method with a point per site (colored by ecosystem,
   panels for both OpenET and flux tower), example monthly timeseries for
   three representative sites (best/median/worst, ranked by calibrated-Oudin
   RMSE vs. OpenET), an ecosystem RMSE comparison, and a spatial map of which
   variant matches OpenET best at each site. Now also evaluates a fourth
   variant (per-ecosystem-calibrated Oudin, from notebook 04) and restricts
   every calibrated-vs-OpenET comparison to a genuine 2022–2023 held-out test
   period, reporting mean/median/worst-decile robustness across sites rather
   than the mean alone. Read-only (no WBM run or GEE access) — just reads
   notebooks 02/03/04/06's cached outputs. See
   `reports/pet_comparison_report.html` for a narrative summary of results.
8. (`notebook 08`) A side analysis, not part of the main pipeline: reuses the
   exact same cached results as notebook 07, filtered down to just the sites
   within 50 km of a park. Since that subset is small (~35 sites), it keeps
   the original full per-site faceted figures (one panel per site) instead
   of notebook 07's condensed versions. No new GEE calls or WBM runs —
   read-only, like notebook 07. **Not yet updated** to notebook 07's held-out,
   four-variant design — still reports the original in-sample, per-site-only
   calibrated comparison for the park-proximity subset.
9. (`notebook 09`) Water-balance implications: since AET is the WBM's primary
   outflow competing with runoff, quantifies how much switching from
   uncalibrated Oudin to per-site-calibrated Oudin or Penman-Monteith would
   change the model's simulated annual runoff, nationwide and by ecosystem.
   Read-only — reuses notebooks 03/04/06's cached daily WBM outputs
   (`RUNOFF` is already a per-day output column in each).

## Repo structure

```
WBM_ET_Exploration/
├── environment.yml           conda environment (see Setup below)
├── wbm/                       NPS WBM model, as an importable Python package
│   ├── met_pet.py              meteorological helpers, daylength, PET (Oudin/Hamon/Penman-Monteith), wind-height correction, GDD/deficit
│   ├── snow_soil.py             precip partitioning, snow/melt, soil moisture, storage reservoir
│   ├── model.py                 main driver (nps_wbm), multi-point wrapper, run_pipeline
│   ├── raster_io.py              DEM/soil/Jennings-temp raster loading + site-parameter extraction
│   └── cache_utils.py            incremental per-site caching (load_and_filter_missing, merge_and_save_cache)
├── scripts/
│   ├── setup_environment.py    run once per session: dependency check, wbm import, load rasters, sanity check
│   └── smoke_test_wbm.py       standalone synthetic-data smoke test (no rasters/data needed)
├── notebooks/                 the analysis pipeline, split into task-based notebooks (run in order)
│   ├── 01_select_parks_towers.ipynb   nationwide site list + park-proximity subset (see below)
│   ├── 02_load_flux_openet_gridmet.ipynb
│   ├── 03_run_wbm.ipynb
│   ├── 04_validate_and_calibrate.ipynb
│   ├── 05_figures_and_export.ipynb
│   ├── 06_oudin_vs_penman_monteith_pet.ipynb   Oudin vs. Penman-Monteith PET/AET comparison (see below)
│   ├── 07_pet_methods_vs_openet.ipynb   both PET methods' WBM AET vs. OpenET, condensed for nationwide scale (see below)
│   ├── 08_park_proximity_side_analysis.ipynb   same comparison, filtered to park-adjacent sites, full per-site facets (see below)
│   ├── 09_runoff_implications.ipynb   water-balance implications of AET method choice (see below)
│   └── exploratory/
│       └── alt_gridmet_thredds.ipynb   untested scratch: THREDDS-based gridMET fetch + raster QA
├── Data/                      NOT tracked in git (see Data below) — inputs, cached pulls, figures
└── docs/                      NOT tracked in git — original class deliverables (report, slides, scope)
```

Each notebook caches its expensive outputs (GEE pulls, WBM runs, per-site
calibration) to CSVs under `Data/`, so re-running a later notebook doesn't
require repeating an earlier notebook's downloads. Within notebooks 02, 03,
04, and 06, this caching is now **incremental per-site** (via
`wbm.cache_utils`): each of those notebooks checks which sites are already in
its cache and only does the expensive work (GEE download, WBM run,
calibration) for sites that are missing, then merges the new rows in. This
means adding a new site to `flux_towers.csv` and re-running the pipeline
doesn't re-download or re-run anything for sites already cached — it only
processes what's new. (This resumes *across* runs, not *within* one: if a
run is interrupted partway through the new sites, that batch's progress is
lost and redone next time — see `wbm/cache_utils.py`'s docstring.)

## Setup

**Enabling the results report link above:** GitHub doesn't render `.html`
files inline, so the report needs [GitHub Pages](https://pages.github.com/)
turned on once: repo **Settings → Pages → Build and deployment → Source:
Deploy from a branch → Branch: `main`, folder: `/ (root)` → Save**. After a
minute it's live at `https://<username>.github.io/<repo>/reports/pet_comparison_report.html`.
It also helps to set the repo's **Website** field (gear icon next to "About"
on the repo's main page) to that same URL — it then shows as a clickable
link right under the repo name for anyone landing on the page.

```bash
conda env create -f environment.yml
conda activate wbm-et-exploration
```

This installs the geospatial stack (rasterio, geopandas, gdal), Google Earth
Engine (`earthengine-api`, `geemap`), and the scientific/plotting packages
the notebooks need. See the comment block at the bottom of `environment.yml`
for a pip-only alternative if you don't use conda.

You'll also need a Google Earth Engine account (free) — `notebooks/01` and
`02` authenticate interactively via `ee.Authenticate()`.

To sanity-check the WBM model code itself (no data or GEE needed):

```bash
python scripts/smoke_test_wbm.py
```

## Data

`Data/` is git-ignored — it's ~1.5 GB and everything in it is either
downloaded or regenerated by the notebooks themselves:

* **NPS park boundaries** — downloaded on demand from the NPS DataStore API
  (see `get_park_system()` in `notebooks/01_select_parks_towers.ipynb`).
* **Flux tower ET** (`Data/flux_ET_dataset/`) — the public Dryad/AmeriFlux
  "Post-processed CONUS eddy flux ET dataset" (Volk et al., 2023a,
  [10.5281/zenodo.7636781](https://doi.org/10.5281/zenodo.7636781));
  download separately and place it at `Data/flux_ET_dataset/` (see its own
  `README.md` once downloaded for the file layout).
* **OpenET / gridMET** — pulled via Google Earth Engine in
  `notebooks/02_load_flux_openet_gridmet.ipynb`, cached to
  `Data/open_et/`, `Data/open_et_gridcell/`, and `Data/gridmet_cache/`.
* **WBM rasters** (`Data/wbm_rasters/`) — elevation, soil water storage
  capacity, and Jennings temperature climatology GeoTIFFs. Used by
  `wbm.raster_io.load_wbm_rasters()`/`extract_point_params()` to sample
  per-site model parameters. These cover the full CONUS extent (confirmed
  against the nationwide flux tower list), not just the original six parks,
  so no new rasters are needed for the nationwide expansion.
* **Notebook 07 vs. 08 figure outputs** — notebook 07 (nationwide, condensed
  figures) saves to `Data/pet_comparison/`; notebook 08 (park-proximity,
  full per-site facets) saves to `Data/pet_comparison_near_parks/` so the
  two don't overwrite each other.

## Known limitations

* Every notebook (01→02→03→04→05, and 06→07/08) reloads its prerequisites
  from `Data/` at the top of the notebook, so each can be opened in a fresh
  kernel as long as the earlier ones have been run at least once. (Notebook
  04's calibration outputs — `OBJECTIVE`, `eco_metrics`, and the
  `metrics_*_vs_*` tables notebook 05 needs — are cached to
  `Data/gridmet_cache/` by a cell at the end of notebook 04.)
* `reports/pet_comparison_report.html` reflects the nationwide, held-out,
  four-variant comparison (notebooks 04/07); notebook 08's park-proximity
  side analysis has not yet been brought up to that same design (see item 8
  above) and notebook 09's runoff-implications findings aren't folded into
  the main report yet.
* The incremental per-site caching (`wbm.cache_utils`, used in notebooks 02,
  03, 04, and 06) resumes *across* runs but doesn't checkpoint *within* one:
  if a run is interrupted partway through the new/missing sites, that
  batch's progress is lost and gets redone on the next run.

## Notes on reproducibility

* Two duplicate 120 MB park-boundary GeoJSON exports from the original
  class project were consolidated into a single
  `Data/geo_data/nps_unit_boundaries.geojson` (not tracked in git, per
  above) — both were byte-identical in content.
* `docs/` holds the original report, slides, project scope, and a
  reference PDF from the class project. It's excluded from git since the
  reference PDF is a copyrighted journal article; the analysis and results
  are documented in the notebooks themselves.
