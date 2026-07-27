# %% [markdown]
# # NPS Water Balance Model — Session Setup
#
# Run the cells in this file **once at the top of every WBM session** (or
# `%run scripts/setup_environment.py` from a notebook under `notebooks/`)
# before any analysis cells. They will:
# 1. Verify all required packages are installed (see `environment.yml`).
# 2. Import the WBM model from the `wbm` package.
# 3. Load and cache the WBM rasters from `Data/wbm_rasters/`.
# 4. Print a quick sanity-check confirming everything is ready.
#
# This replaces the old class-project `nps_wbm_setup.py`, updated to import
# from the `wbm` package (split out of the original single-file
# `nps_wbm_functions.py`) instead of a hard-coded class-assignment path.

# %% [markdown]
# ## Cell 1 — Dependency check
#
# Checks that the core packages are importable. Geospatial/GEE packages
# (rasterio, geopandas, earthengine-api, geemap, ...) are managed via
# `environment.yml` — see the repo README for setup instructions.

# %%
import importlib
import importlib.util
import sys

REQUIRED = {
    "numpy":    "numpy",
    "pandas":   "pandas",
    "scipy":    "scipy",
    "rasterio": "rasterio",
}

missing = [pip_name for mod_name, pip_name in REQUIRED.items()
           if importlib.util.find_spec(mod_name) is None]

if missing:
    print(f"WARNING: missing packages: {missing}")
    print("Create/update the conda environment first:\n")
    print("  conda env create -f environment.yml")
    print("  conda activate wbm-et-exploration\n")
else:
    print("All required packages are available.")
    print(f"  Python  : {sys.version.split()[0]}")
    for mod_name in REQUIRED:
        mod = importlib.import_module(mod_name)
        ver = getattr(mod, "__version__", "unknown")
        print(f"  {mod_name:<10}: {ver}")

# %% [markdown]
# ## Cell 2 — Import the WBM model
#
# The `wbm` package lives at the repo root (`wbm/`). If you're running this
# from a notebook under `notebooks/`, make sure the repo root is on
# `sys.path` first (the notebooks in this repo already do this).

# %%
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent if "__file__" in dir() else Path.cwd()
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import wbm
importlib.reload(wbm)

from wbm import (
    # ── Raster utilities ──────────────────────────────────────────────────
    load_wbm_rasters,       # load & cache DEM / soil / Jennings rasters
    extract_point_params,   # sample site params at point locations
    run_pipeline,           # one-call convenience: load → extract → run

    # ── Core model ───────────────────────────────────────────────────────
    nps_wbm,                # single-point daily water balance driver
    run_nps_wbm_points,     # multi-point / multi-GCM wrapper

    # ── Component functions (available if needed for custom workflows) ────
    get_freeze,             # rain/snow partitioning factor
    get_rain, get_snow,     # rainfall and snowfall
    get_melt,               # Hock degree-day snowmelt
    get_snowpack,           # snowpack accumulation
    get_ablation,           # snow sublimation / vapor loss
    get_soil,               # soil water content
    get_d_soil,             # daily change in SWC
    get_aet,                # actual evapotranspiration
    get_storage,            # linear storage reservoir (CSU addition)
    get_oudin_pet,          # Oudin PET with topographic heat-load
    get_hamon_pet,          # Hamon PET
    get_penman_monteith_pet,# FAO-56 Penman-Monteith PET
    get_daylength,          # astronomical daylength (hours)
    get_gdd,                # growing degree days
    get_deficit,            # climatic water deficit (PET − AET)

    # ── Low-level helpers (rarely called directly) ────────────────────────
    get_svp, actual_vp, atm_press, psyc_constant,
    vapor_curve, clear_sky_rad, outgoing_rad,
)

print(f"WBM functions loaded from: {REPO_ROOT / 'wbm'}")

# %% [markdown]
# ## Cell 3 — Load WBM rasters
#
# Reads the three required GeoTIFFs from `Data/wbm_rasters/`, derives slope
# and aspect from the DEM (Horn's 8-neighbour method, matching
# `terra::terrain(neighbors = 8)`), and caches everything in memory.
#
# **Only needs to run once per session.** After this, `extract_point_params()`
# works without any path arguments.
#
# | File | Content | Units |
# |---|---|---|
# | `elevation_cropped.tif` | DEM | metres |
# | `water_storage.tif` | Soil water storage capacity | cm (auto-converted → mm) |
# | `merged_jennings2.tif` | Jennings temperature climatology | °C |

# %%
# ── Path to raster directory ──────────────────────────────────────────────────
# Relative to the repo root. If calling from a notebook under notebooks/,
# use "../Data/wbm_rasters" instead.
RASTER_DIR = str(REPO_ROOT / "Data" / "wbm_rasters")

# ── Target CRS ────────────────────────────────────────────────────────────────
# Set to None to keep the native CRS of the rasters (EPSG:4326 for NPS data).
# Set to an EPSG string (e.g. "EPSG:26913") to reproject all rasters before
# sampling — useful when your point coordinates are in a projected CRS.
TARGET_CRS = None   # e.g. "EPSG:4326"

try:
    load_wbm_rasters(raster_dir=RASTER_DIR, target_crs=TARGET_CRS)
except ImportError as e:
    print(f"WARNING: rasters not loaded: {e}")
    print("   Install rasterio/scipy and re-run this cell before calling extract_point_params().")
except FileNotFoundError as e:
    print(f"WARNING: raster file not found:\n   {e}")
    print(f"   Check that RASTER_DIR = '{RASTER_DIR}' is correct and that Data/wbm_rasters/ exists locally")
    print("   (Data/ is git-ignored — see README for how to regenerate/download it).")

# %% [markdown]
# ## Cell 4 — Global model settings
#
# Set your default WBM run parameters here. These are passed through to
# `nps_wbm()` / `run_nps_wbm_points()` / `run_pipeline()` in later cells.
# Override any of them at call-time as needed.

# %%
import numpy as np
import pandas as pd

# ── PET method ────────────────────────────────────────────────────────────────
# One of: "Oudin"  (default, temperature-based, topographic heat-load adjusted)
#         "Hamon"  (temperature + daylength)
#         "Penman-Monteith"  (requires tmax, tmin, and ideally RH / wind data)
PET_METHOD = "Oudin"

# ── Snowmelt ──────────────────────────────────────────────────────────────────
# Hock (2003) degree-day melt factor (mm °C⁻¹ day⁻¹).
# Hock reports ~2.5 for Gooseberry Creek, UT; NPS default is 4.
HOCK_COEF = 4.0

# ── Initial conditions ────────────────────────────────────────────────────────
SNOWPACK_INIT = 0.0   # mm SWE
SOIL_INIT     = 0.0   # mm

# ── CSU additions ─────────────────────────────────────────────────────────────
DIRECT_FRAC  = 0.0    # fraction of rainfall routed directly to runoff (0–1)
RETURN_RATE  = 1.0    # fraction of storage reservoir released per day (0–1]
PET_MULT     = 1.0    # multiplicative PET bias correction
SOIL_MULT    = 1.0    # multiplicative adjustment to SWC_Max

# ── Misc ──────────────────────────────────────────────────────────────────────
SHADE_COEFF  = 1.0    # canopy shading coefficient for Oudin PET (0–1)
T_BASE       = 0.0    # base temperature for growing degree days (°C)
TO_INCHES    = True   # True → output fluxes in inches; False → mm

print("Model settings configured:")
print(f"  PET method   : {PET_METHOD}")
print(f"  Hock coef    : {HOCK_COEF} mm °C⁻¹ day⁻¹")
print(f"  Direct frac  : {DIRECT_FRAC}")
print(f"  Return rate  : {RETURN_RATE}")
print(f"  PET mult     : {PET_MULT}   |  Soil mult: {SOIL_MULT}")
print(f"  Output units : {'inches' if TO_INCHES else 'mm'}")

# %% [markdown]
# ## Cell 5 — Verify setup (quick sanity check)
#
# Runs the model for one synthetic year at a single point to confirm the full
# stack (imports → raster cache → model) is working before you connect real
# climate data. (See also `scripts/smoke_test_wbm.py` for a standalone
# version of this check.)

# %%
rng = np.random.default_rng(0)
_n  = 365
_dates = pd.date_range("2000-01-01", periods=_n)
_doy   = _dates.dayofyear.to_numpy(float)

_test_climate = pd.DataFrame({
    "date":    _dates,
    "x":       -105.5,          # lon — update to match your study area
    "y":        40.0,           # lat
    "ppt_mm":  rng.exponential(3.0, _n),
    "tmean_C": 8 * np.sin(2 * np.pi * (_doy - 80) / 365) + 5 + rng.normal(0, 2, _n),
    "GCM":     "sanity_check",
})

# Use hard-coded params so the test doesn't depend on the rasters being loaded
_test_params = {"Elev": 2400, "Slope": 10, "Aspect": 180,
                "SWC_Max": 150, "J_Temp": 1.5}

_test_result = nps_wbm(
    _test_climate, _test_params,
    pet_method    = PET_METHOD,
    hock_coef     = HOCK_COEF,
    direct_frac   = DIRECT_FRAC,
    return_rate   = RETURN_RATE,
    pet_mult      = PET_MULT,
    soil_mult     = SOIL_MULT,
    shade_coeff   = SHADE_COEFF,
    t_base        = T_BASE,
    to_inches     = False,          # mm for the sanity check
)

_unit = "mm"
print("Sanity check passed — annual water balance totals:")
print(f"  {'Variable':<14}  {'Annual total':>14}")
print(f"  {'-'*30}")
for _col in ["ppt_mm", "RAIN", "SNOW", "MELT", "AET", "RUNOFF", "D"]:
    print(f"  {_col:<14}  {_test_result[_col].sum():>12.1f} {_unit}")

print(f"\n  Peak snowpack : {_test_result['PACK'].max():.1f} {_unit}")
print(f"  Max soil SWC  : {_test_result['SOIL'].max():.1f} {_unit}")
print(f"\nSetup complete — ready to run NPS WBM.\n")

# %% [markdown]
# ---
# ## Quick-reference: key function signatures
#
# ```python
# # ── Extract site params at your points from the loaded rasters ────────────
# point_params_df = extract_point_params(points_df)
# # points_df needs columns: x (lon), y (lat)
# # Returns:  Elev, Slope, Aspect, SWC_Max, J_Temp  added to points_df
#
# # ── Run for a single point ────────────────────────────────────────────────
# result = nps_wbm(
#     daily_df     = climate_df,        # date, x, y, ppt_mm, tmean_C [, GCM]
#     point_params = point_params_df.iloc[0].to_dict(),
#     pet_method   = PET_METHOD,
#     **{k: v for k, v in globals().items()
#        if k in ("hock_coef","direct_frac","return_rate","pet_mult",
#                 "soil_mult","shade_coeff","t_base","to_inches")},
# )
#
# # ── Run for multiple points / GCMs ───────────────────────────────────────
# results = run_nps_wbm_points(
#     climate_data    = climate_df,
#     point_params_df = point_params_df,
#     pet_method      = PET_METHOD,
#     aggregate       = True,           # False → keep per-point rows
#     ret_final_cond  = False,          # True → chain into next time period
# )
#
# # ── One-call pipeline (load → extract → run) ─────────────────────────────
# results = run_pipeline(
#     climate_data = climate_df,
#     points_df    = points_df,
#     raster_dir   = RASTER_DIR,
#     pet_method   = PET_METHOD,
# )
# ```
