"""
NPS Water Balance Model - main driver, multi-point wrapper, and end-to-end
pipeline.

Ported from the R implementation originally developed by ARC (10/30/2019)
with updates by MGS (4/14/2021) and CSU additions. Based on Dave Thoma's
water balance Excel spreadsheet model.

Typical usage
-------------
    import pandas as pd
    from wbm import nps_wbm, run_nps_wbm_points

    # --- single point ---
    results = nps_wbm(
        daily_df   = climate_df,        # DataFrame with date, ppt_mm, tmean_C, x, y
        point_params = {
            "Elev":    2400,            # m
            "Slope":   15,              # degrees
            "Aspect":  180,             # degrees
            "SWC_Max": 150,             # mm
            "J_Temp":  1.5,             # °C  (Jennings climatology)
        },
        PET_Method = "Oudin",
    )

    # --- multiple points / GCMs ---
    results = run_nps_wbm_points(
        climate_data = climate_df,      # must include 'x', 'y', 'GCM' columns
        point_params_df = params_df,    # one row per point
        PET_Method = "Oudin",
    )

Parameter extraction from rasters
----------------------------------
    If you have GeoTIFF rasters for elevation, soil, slope, aspect, and
    Jennings temperature, use ``wbm.raster_io.extract_point_params()``
    (or the ``run_pipeline()`` convenience wrapper below), which depends
    on ``rasterio`` (optional).
"""

from __future__ import annotations

import warnings
from typing import Optional, Union

import numpy as np
import pandas as pd

from .met_pet import (
    get_daylength,
    get_oudin_pet,
    get_hamon_pet,
    get_penman_monteith_pet,
    get_gdd,
    get_deficit,
)
from .snow_soil import (
    get_freeze,
    get_rain,
    get_snow,
    get_melt,
    get_snowpack,
    get_soil,
    get_d_soil,
    get_aet,
    get_storage,
)
from .raster_io import extract_point_params, _DEFAULT_RASTER_DIR

# ──────────────────────────────────────────────────────────────────────────────
# 9.  MAIN MODEL DRIVER
# ──────────────────────────────────────────────────────────────────────────────

def nps_wbm(daily_df: pd.DataFrame,
            point_params: dict,
            direct_frac: float = 0.0,
            return_rate: float = 1.0,
            pet_mult: float = 1.0,
            soil_mult: float = 1.0,
            pet_method: str = "Oudin",
            hock_coef: float = 4.0,
            snowpack_init: float = 0.0,
            soil_init: float = 0.0,
            shade_coeff: float = 1.0,
            t_base: float = 0.0,
            to_inches: bool = True) -> pd.DataFrame:
    """
    NPS Water Balance Model – single-point daily driver.

    Parameters
    ----------
    daily_df : pd.DataFrame
        Daily climate records. Required columns:

        * ``date``    – datetime or date-like
        * ``ppt_mm``  – precipitation (mm day⁻¹)
        * ``tmean_C`` – mean air temperature (°C)
        * ``x``       – longitude (decimal degrees)
        * ``y``       – latitude  (decimal degrees)

        Optional (needed for Penman-Monteith):
        ``tmax_C``, ``tmin_C``, ``RHmax``, ``RHmin``, ``vp``

        Any additional columns (e.g., ``GCM``) are preserved in the output.

    point_params : dict
        Site parameters with keys:

        * ``Elev``    – elevation (m)
        * ``Slope``   – slope (degrees)
        * ``Aspect``  – aspect (degrees, 0–360)
        * ``SWC_Max`` – max soil water capacity (mm)
        * ``J_Temp``  – Jennings temperature climatology (°C)

        The latitude is taken from ``daily_df["y"]`` so that spatially
        varying lat is supported automatically.

    direct_frac : float
        Fraction of rainfall routed directly to runoff (0–1). Default 0.
    return_rate : float
        Daily fraction of storage reservoir released to runoff (0–1].
        Default 1 (full same-day release).
    pet_mult : float
        Multiplicative PET bias adjustment. Default 1.
    soil_mult : float
        Multiplier on SWC_Max. Default 1.
    pet_method : str
        One of ``"Oudin"``, ``"Hamon"``, or ``"Penman-Monteith"``.
    hock_coef : float
        Degree-day melt factor (mm °C⁻¹ day⁻¹). Default 4.
    snowpack_init : float
        Initial snowpack (mm SWE). Default 0.
    soil_init : float
        Initial soil water content (mm). Default 0.
    shade_coeff : float
        Shade coefficient for Oudin PET (0–1). Default 1.
    t_base : float
        Base temperature for GDD (°C). Default 0.
    to_inches : bool
        If True (default), convert all flux outputs from mm to inches
        (temperature columns remain °C).

    Returns
    -------
    pd.DataFrame
        Wide-format daily water-balance results. Columns:

        ``date, x, y, [GCM], ppt_mm,``
        ``RAIN, SNOW, MELT, PACK, W, PET, PET_mod,``
        ``SOIL, DSOIL, AET, D (deficit),``
        ``STORAGE_ADD, STORAGE_RELEASE, STORAGE_REMAIN,``
        ``RUNOFF, GDD, tmean_C``

        Flux columns are in inches when ``to_inches=True``, mm otherwise.
    """
    df = daily_df.copy()
    df["date"] = pd.to_datetime(df["date"])
    df = df.sort_values("date").reset_index(drop=True)

    lat = df["y"].to_numpy(dtype=float)
    elev = float(point_params["Elev"])
    slope = float(point_params["Slope"])
    aspect = float(point_params["Aspect"])
    swc_max = float(point_params["SWC_Max"]) * soil_mult
    j_temp = float(point_params["J_Temp"])

    ppt = df["ppt_mm"].to_numpy(dtype=float)
    tmean = df["tmean_C"].to_numpy(dtype=float)
    doy = df["date"].dt.dayofyear.to_numpy(dtype=float)
    jtemp_arr = np.full(len(df), j_temp)

    # ── Day length ──────────────────────────────────────────────────────────
    df["daylength"] = get_daylength(df["date"], lat)

    # ── Rain / snow partitioning ────────────────────────────────────────────
    freeze = get_freeze(jtemp_arr, tmean)
    rain = get_rain(ppt, freeze)
    snow = get_snow(ppt, freeze)

    # ── Snowmelt & snowpack ─────────────────────────────────────────────────
    melt = get_melt(tmean, jtemp_arr, hock_coef, snow, sp0=snowpack_init)
    pack = get_snowpack(jtemp_arr, snow, melt, sp0=snowpack_init)

    # ── Direct runoff fraction (CSU addition) ───────────────────────────────
    direct = rain * direct_frac
    w = melt + rain - direct              # water available to soil

    # ── PET ─────────────────────────────────────────────────────────────────
    pet_method = pet_method.strip()
    if pet_method == "Hamon":
        pet = get_hamon_pet(tmean, df["daylength"].to_numpy())
    elif pet_method == "Penman-Monteith":
        tmax = df.get("tmax_C", pd.Series(tmean + 5)).to_numpy(dtype=float)
        tmin = df.get("tmin_C", pd.Series(tmean - 5)).to_numpy(dtype=float)
        rhmax = df["RHmax"].to_numpy() if "RHmax" in df else None
        rhmin = df["RHmin"].to_numpy() if "RHmin" in df else None
        vp_col = df["vp"].to_numpy() if "vp" in df else None
        # Optional daily wind speed at 2 m (e.g. gridMET 'vs' band corrected
        # via get_wind_speed_2m()). Falls back to the function's constant
        # default (2.2 m/s) if no 'wind_ms' column is supplied — previously
        # this branch never passed wind through at all, so any wind_ms
        # column was silently ignored.
        #
        # NOTE: a present-but-NaN 'wind_ms' column (e.g. a site cached before
        # wind was added to a climate pull, then never re-fetched — this
        # happened for ~35/81 sites in this repo's own gridMET cache) used to
        # pass NaN straight into get_penman_monteith_pet(), which silently
        # produced NaN PET/AET for the entire record. Downstream monthly
        # aggregation (pandas .sum() on an all-NaN group) then silently
        # reported that as AET = 0.0 mm/month with no error or warning. The
        # check below falls back to the 2.2 m/s default on a per-row basis
        # wherever wind is actually missing, and prints a one-line warning so
        # a data gap like that can't disappear silently again.
        if "wind_ms" in df:
            wind_arg = df["wind_ms"].to_numpy(dtype=float)
            n_missing = np.isnan(wind_arg).sum()
            if n_missing > 0:
                site_label = df["site"].iloc[0] if "site" in df else "unknown site"
                print(f"WARNING: {n_missing}/{len(wind_arg)} 'wind_ms' values are NaN "
                      f"for {site_label} — falling back to the 2.2 m/s default for "
                      "those rows instead of propagating NaN into PET/AET.")
                wind_arg = np.where(np.isnan(wind_arg), 2.2, wind_arg)
        else:
            wind_arg = 2.2
        pet = get_penman_monteith_pet(tmax, tmin, doy, elev,
                                      lat[0], rhmax, rhmin, vp_col,
                                      wind=wind_arg)
    elif pet_method == "Oudin":
        pet = get_oudin_pet(doy, lat[0], pack, tmean, slope, aspect, shade_coeff)
    else:
        raise ValueError(f"Unknown PET method: '{pet_method}'. "
                         "Choose 'Oudin', 'Hamon', or 'Penman-Monteith'.")

    pet_mod = pet * pet_mult

    # ── Soil moisture ────────────────────────────────────────────────────────
    soil = get_soil(w, pet_mod, swc_max, swc0=soil_init)
    d_soil = get_d_soil(soil, swc0=soil_init)
    aet = get_aet(w, pet_mod, soil, swc0=soil_init)

    # ── Storage reservoir ────────────────────────────────────────────────────
    storage_add = w - aet - d_soil
    stor = get_storage(storage_add, return_rate)

    # ── Runoff & deficit ─────────────────────────────────────────────────────
    runoff = stor["storage_release"] + direct
    deficit = get_deficit(pet_mod, aet)

    # ── Growing degree days ──────────────────────────────────────────────────
    gdd = get_gdd(tmean, t_base)

    # ── Assemble output ───────────────────────────────────────────────────────
    out = df[["date", "x", "y"]].copy()
    if "GCM" in df.columns:
        out["GCM"] = df["GCM"]

    flux_cols = {
        "ppt_mm": ppt,
        "RAIN": rain,
        "SNOW": snow,
        "MELT": melt,
        "PACK": pack,
        "W": w,
        "PET": pet,
        "PET_mod": pet_mod,
        "SOIL": soil,
        "DSOIL": d_soil,
        "AET": aet,
        "D": deficit,
        "STORAGE_ADD": storage_add,
        "STORAGE_RELEASE": stor["storage_release"],
        "STORAGE_REMAIN": stor["storage_remain"],
        "RUNOFF": runoff,
    }

    scale = 1.0 / 25.4 if to_inches else 1.0
    for col, arr in flux_cols.items():
        out[col] = arr * scale

    # temperature stays in °C
    out["GDD"] = gdd
    out["tmean_C"] = tmean

    return out


# ──────────────────────────────────────────────────────────────────────────────
# 10. MULTI-POINT / MULTI-GCM WRAPPER
# ──────────────────────────────────────────────────────────────────────────────

def run_nps_wbm_points(climate_data: pd.DataFrame,
                        point_params_df: pd.DataFrame,
                        direct_frac: float = 0.0,
                        return_rate: float = 1.0,
                        pet_mult: float = 1.0,
                        soil_mult: float = 1.0,
                        pet_method: str = "Oudin",
                        hock_coef: float = 4.0,
                        snowpack_init: float = 0.0,
                        soil_init: float = 0.0,
                        shade_coeff: float = 1.0,
                        t_base: float = 0.0,
                        to_inches: bool = True,
                        aggregate: bool = True,
                        ret_final_cond: bool = False,
                        use_final: bool = False,
                        area_weights: Optional[np.ndarray] = None,
                        ) -> Union[pd.DataFrame, tuple[pd.DataFrame, pd.DataFrame]]:
    """
    Multi-point / multi-GCM Water Balance Model wrapper.

    Loops over all (point, GCM) combinations, runs ``nps_wbm()``, and
    optionally aggregates results spatially.

    Parameters
    ----------
    climate_data : pd.DataFrame
        Daily climate for all points and GCMs. Required columns:

        * ``x``      – longitude (decimal degrees)
        * ``y``      – latitude  (decimal degrees)
        * ``date``   – datetime-like
        * ``ppt_mm`` – precipitation (mm day⁻¹)
        * ``tmean_C``– mean temperature (°C)
        * ``GCM``    – model/scenario identifier (str)

    point_params_df : pd.DataFrame
        One row per point with columns matching ``nps_wbm()``'s
        ``point_params`` dict keys plus ``x`` and ``y`` for merging:

        ``x, y, Elev, Slope, Aspect, SWC_Max, J_Temp``

        When ``use_final=True``, must also contain ``final_soil`` and
        ``final_snowpack`` columns (populated by a prior run with
        ``ret_final_cond=True``).

    direct_frac, return_rate, pet_mult, soil_mult,
    pet_method, hock_coef, snowpack_init, soil_init,
    shade_coeff, t_base, to_inches :
        Passed through to ``nps_wbm()`` – see its docstring.

    aggregate : bool
        If True (default), return daily means across all points (or
        area-weighted means when ``area_weights`` is supplied).
        If False, return per-point rows with ``x`` and ``y`` retained.

    ret_final_cond : bool
        If True, also return the updated ``point_params_df`` with columns
        ``final_soil`` and ``final_snowpack`` filled in from the last
        time step of each point's run.

    use_final : bool
        If True, initialise each point from ``point_params_df``'s
        ``final_soil`` / ``final_snowpack`` columns (requires those
        columns to exist).

    area_weights : array-like, optional
        One weight per row of ``point_params_df``. When supplied with
        ``aggregate=True``, weighted means are used. Weights are
        normalised internally so they need not sum to 1.

    Returns
    -------
    pd.DataFrame
        Daily water-balance results (wide format, columns as in
        ``nps_wbm()`` output). When ``aggregate=True`` the spatial
        dimension is collapsed; otherwise each row includes ``x`` and ``y``.

    tuple[pd.DataFrame, pd.DataFrame]
        Only when ``ret_final_cond=True``: ``(results_df, point_params_df)``
        where ``point_params_df`` now contains final soil/snowpack states.
    """
    params = point_params_df.copy().reset_index(drop=True)

    # Initialise final condition columns if needed
    if ret_final_cond and "final_soil" not in params.columns:
        params["final_soil"] = 0.0
    if ret_final_cond and "final_snowpack" not in params.columns:
        params["final_snowpack"] = 0.0

    all_gcms = climate_data["GCM"].unique()
    gcm_results = []

    print("Running NPS WBM...")
    for gcm in all_gcms:
        print(f"  GCM: {gcm}")
        point_results = []

        for idx, row in params.iterrows():
            print(f"    Point {idx + 1}/{len(params)}")
            px, py = float(row["x"]), float(row["y"])

            # Per-point init overrides when chaining runs
            sp_init = float(row["final_snowpack"]) if use_final else snowpack_init
            s_init = float(row["final_soil"]) if use_final else soil_init

            # Subset climate to this point and GCM
            mask = (
                np.isclose(climate_data["x"].to_numpy(dtype=float), px)
                & np.isclose(climate_data["y"].to_numpy(dtype=float), py)
                & (climate_data["GCM"] == gcm)
            )
            pt_climate = climate_data.loc[mask].copy()

            if pt_climate.empty:
                warnings.warn(
                    f"No climate data for point ({px}, {py}) / GCM '{gcm}'. Skipping.",
                    stacklevel=2,
                )
                continue

            pp = row.to_dict()
            result = nps_wbm(
                daily_df=pt_climate,
                point_params=pp,
                direct_frac=direct_frac,
                return_rate=return_rate,
                pet_mult=pet_mult,
                soil_mult=soil_mult,
                pet_method=pet_method,
                hock_coef=hock_coef,
                snowpack_init=sp_init,
                soil_init=s_init,
                shade_coeff=shade_coeff,
                t_base=t_base,
                to_inches=to_inches,
            )
            result["GCM"] = gcm

            # Save final conditions
            if ret_final_cond:
                scale = 25.4 if to_inches else 1.0   # back to mm
                params.at[idx, "final_soil"] = float(result["SOIL"].iloc[-1]) * scale
                params.at[idx, "final_snowpack"] = float(result["PACK"].iloc[-1]) * scale

            # Attach area weight if provided
            if area_weights is not None:
                result["area_weight"] = float(
                    np.asarray(area_weights, dtype=float)[idx]
                )

            point_results.append(result)

        if not point_results:
            continue

        combined = pd.concat(point_results, ignore_index=True)

        if aggregate:
            flux_cols = [c for c in combined.columns
                         if c not in ("date", "x", "y", "GCM", "area_weight")]
            grp = combined.groupby(["date", "GCM"])

            if area_weights is not None:
                # Normalised weighted mean per day
                def w_mean(grp_df):
                    w = grp_df["area_weight"]
                    w_norm = w / w.sum()
                    return (grp_df[flux_cols].multiply(w_norm.values, axis=0)
                            .sum())
                agg = grp.apply(w_mean).reset_index()
            else:
                agg = grp[flux_cols].mean().reset_index()

            gcm_results.append(agg)
        else:
            gcm_results.append(combined)

    if not gcm_results:
        raise RuntimeError("No results were produced. Check climate_data alignment with point_params_df.")

    final_df = pd.concat(gcm_results, ignore_index=True)

    # Rename to match R output conventions
    rename_map = {
        "RUNOFF": "runoff_wbm",
        "RAIN": "rain_wbm",
        "SNOW": "snow_wbm",
        "MELT": "melt_wbm",
        "AET": "aet_wbm",
        "ppt_mm": "precip_wbm",
        "PET_mod": "pet_wbm",
        "STORAGE_ADD": "excess_water_wbm",
        "SOIL": "soil_wbm",
        "PACK": "snowpack_wbm",
    }
    final_df = final_df.rename(columns=rename_map)

    if ret_final_cond:
        return final_df, params
    return final_df


# ──────────────────────────────────────────────────────────────────────────────
# 12. CONVENIENCE: full pipeline in one call
# ──────────────────────────────────────────────────────────────────────────────

def run_pipeline(climate_data: pd.DataFrame,
                 points_df: pd.DataFrame,
                 raster_dir: str = _DEFAULT_RASTER_DIR,
                 target_crs: Optional[str] = None,
                 **wbm_kwargs) -> pd.DataFrame:
    """
    End-to-end convenience wrapper.

    1. Loads rasters (if not already cached).
    2. Extracts site parameters at each point.
    3. Runs ``run_nps_wbm_points()``.

    Parameters
    ----------
    climate_data : pd.DataFrame
        Daily climate (``x``, ``y``, ``date``, ``ppt_mm``, ``tmean_C``,
        ``GCM``).
    points_df : pd.DataFrame
        Point locations (``x``, ``y``).
    raster_dir : str
        Directory with WBM GeoTIFF files.
    target_crs : str or None
        CRS to reproject rasters to. Default keeps native CRS.
    **wbm_kwargs
        Any keyword argument accepted by ``run_nps_wbm_points()``
        (e.g. ``pet_method``, ``hock_coef``, ``to_inches``).

    Returns
    -------
    pd.DataFrame
        Daily water-balance results.

    Example
    -------
    >>> results = run_pipeline(
    ...     climate_data = my_climate_df,
    ...     points_df    = my_points_df,
    ...     raster_dir   = "Data/wbm_rasters",
    ...     pet_method   = "Oudin",
    ...     to_inches    = False,
    ... )
    """
    point_params_df = extract_point_params(points_df,
                                            raster_dir=raster_dir,
                                            target_crs=target_crs)
    return run_nps_wbm_points(climate_data=climate_data,
                               point_params_df=point_params_df,
                               **wbm_kwargs)
