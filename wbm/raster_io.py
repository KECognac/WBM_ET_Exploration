"""
Raster loading and geospatial parameter extraction for the NPS Water
Balance Model (replaces R's ``load_wbm_rasters`` / ``get_nps_wbm_params``).

Design note
-----------
R used global ``<<-`` assignment so every downstream function could see the
rasters. Here we use a module-level ``_RASTERS`` dict that is populated once
by ``load_wbm_rasters()`` and then read by ``extract_point_params()``. Call
``load_wbm_rasters()`` once at the top of your script; after that,
``extract_point_params()`` works without any path arguments.

Raster files expected (matching original R paths, updated directory):
    Data/wbm_rasters/elevation_cropped.tif   – DEM in metres
    Data/wbm_rasters/water_storage.tif       – soil storage in cm (×10 → mm)
    Data/wbm_rasters/merged_jennings2.tif    – Jennings temp climatology (°C)

Slope & aspect are derived from the DEM using Horn's 8-neighbour weighted
finite-difference method – the same algorithm used by ``terra::terrain``
with ``neighbors = 8`` in the original R code.

This is the only module in the ``wbm`` package that requires ``rasterio``
and ``scipy`` — those imports are deferred (via ``_require_rasterio()`` /
``_require_scipy()``) so the rest of the package (pure WBM math) works even
if they aren't installed.
"""

from __future__ import annotations

import warnings
from typing import Optional

import numpy as np
import pandas as pd

_RASTERS: dict = {}          # module-level cache; populated by load_wbm_rasters()

# Default raster directory (relative to the current working directory when a
# script/notebook is run). Notebooks under notebooks/ should pass
# raster_dir="../Data/wbm_rasters" explicitly since they run one level
# deeper than the repo root.
_DEFAULT_RASTER_DIR = "Data/wbm_rasters"


def _require_rasterio():
    """Import rasterio (and rasterio.warp) or raise a clear error."""
    try:
        import rasterio                          # noqa: F401
        import rasterio.warp                     # noqa: F401
        import rasterio.crs                      # noqa: F401
        return rasterio
    except ImportError:
        raise ImportError(
            "rasterio is required for raster I/O.\n"
            "Install it with:  pip install rasterio"
        )


def _require_scipy():
    """Import scipy.ndimage or raise a clear error."""
    try:
        from scipy.ndimage import convolve       # noqa: F401
        return convolve
    except ImportError:
        raise ImportError(
            "scipy is required for Horn's slope/aspect computation.\n"
            "Install it with:  pip install scipy"
        )


def _reproject_raster(src_array: np.ndarray,
                       src_transform,
                       src_crs,
                       dst_crs,
                       resampling_str: str = "bilinear"):
    """
    Reproject a raster array to a new CRS.

    Returns (reprojected_array, new_transform, new_crs).
    Equivalent to terra::project() in R.
    """
    rasterio = _require_rasterio()
    from rasterio.warp import reproject, Resampling, calculate_default_transform

    resampling_map = {
        "nearest":  Resampling.nearest,
        "bilinear": Resampling.bilinear,
        "cubic":    Resampling.cubic,
    }
    resampling = resampling_map.get(resampling_str, Resampling.bilinear)

    height, width = src_array.shape
    dst_transform, dst_width, dst_height = calculate_default_transform(
        src_crs, dst_crs, width, height,
        # use src_transform to infer extent
        left   = src_transform.c,
        top    = src_transform.f,
        right  = src_transform.c + src_transform.a * width,
        bottom = src_transform.f + src_transform.e * height,
    )
    dst_array = np.empty((dst_height, dst_width), dtype=np.float64)
    reproject(
        source=src_array.astype(np.float64),
        destination=dst_array,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        resampling=resampling,
    )
    return dst_array, dst_transform, dst_crs


def _compute_slope_aspect_horn(dem_arr: np.ndarray,
                                cell_width: float,
                                cell_height: float
                                ) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute slope (degrees) and aspect (degrees, 0–360) using Horn's
    8-neighbour weighted finite-difference method.

    This replicates terra::terrain(dem, v = c("slope","aspect"),
    unit = "degrees", neighbors = 8).

    Horn (1981) kernels:
        dz/dx = [(-1  0  1)    dz/dy = [( 1  2  1)
                 (-2  0  2)  /           ( 0  0  0)
                 (-1  0  1)] / (8*w)     (-1 -2 -1)] / (8*h)

    Slope  = atan(sqrt(dz/dx² + dz/dy²))
    Aspect = atan2(-dz/dy, dz/dx), rotated to compass bearing 0-360
             where 0/360 = North, 90 = East, 180 = South, 270 = West.

    Parameters
    ----------
    dem_arr : np.ndarray  (2-D, float)
        DEM elevation values.
    cell_width : float
        Horizontal cell size in the same units as elevation (metres or
        degrees – for geographic CRS this should be converted to metres
        before calling if accurate slope is needed).
    cell_height : float
        Vertical cell size (positive value; sign is handled internally).

    Returns
    -------
    slope  : np.ndarray  (degrees, 0–90)
    aspect : np.ndarray  (degrees, 0–360)
    """
    convolve = _require_scipy()

    # Horn's weighted kernels (normalised by 8 so the convolution already
    # accounts for the 8-neighbour weighting)
    kernel_x = np.array([[-1,  0,  1],
                          [-2,  0,  2],
                          [-1,  0,  1]], dtype=float) / 8.0

    kernel_y = np.array([[ 1,  2,  1],
                          [ 0,  0,  0],
                          [-1, -2, -1]], dtype=float) / 8.0

    dz_dx = convolve(dem_arr, kernel_x, mode="nearest") / cell_width
    dz_dy = convolve(dem_arr, kernel_y, mode="nearest") / cell_height

    slope_rad = np.arctan(np.sqrt(dz_dx ** 2 + dz_dy ** 2))
    slope_deg = np.degrees(slope_rad)

    # Aspect: compass bearing where 0 = North.
    # terra convention: atan2(-dz/dy, dz/dx) gives mathematical angle east-of-
    # north; we rotate to compass bearing by: aspect = 90 - angle, then mod 360.
    math_angle = np.degrees(np.arctan2(-dz_dy, dz_dx))
    aspect_deg = (90.0 - math_angle) % 360.0

    return slope_deg, aspect_deg


def load_wbm_rasters(target_crs: Optional[str] = None,
                      raster_dir: str = _DEFAULT_RASTER_DIR) -> None:
    """
    Load and pre-process all WBM raster inputs into the module cache.

    Call this **once** at the start of your script. After calling it,
    use ``extract_point_params()`` without any path arguments.

    Equivalent to R's ``load_wbm_rasters(crs)`` which assigns rasters to
    the global environment with ``<<-``.

    Parameters
    ----------
    target_crs : str or None
        EPSG string (e.g. ``"EPSG:4326"``) or any CRS string accepted by
        rasterio. All rasters are reprojected to this CRS if they differ
        from it. If ``None`` (default), rasters are kept in their native
        CRS (typically EPSG:4326 for the NPS datasets) and no reprojection
        is performed.
    raster_dir : str
        Directory containing the three GeoTIFF files. Defaults to
        ``"Data/wbm_rasters"``; pass ``"../Data/wbm_rasters"`` when calling
        from a notebook under ``notebooks/``.

    Raster files loaded
    -------------------
    ``elevation_cropped.tif``  – DEM (metres)
    ``water_storage.tif``      – Soil water storage (cm; multiplied ×10 → mm)
    ``merged_jennings2.tif``   – Jennings temperature climatology (°C)

    Derived layers
    --------------
    ``slope``  – computed from DEM via Horn's 8-neighbour method (degrees)
    ``aspect`` – computed from DEM via Horn's 8-neighbour method (0–360°)

    Raises
    ------
    FileNotFoundError
        If any of the three required GeoTIFF files cannot be found.
    ImportError
        If ``rasterio`` or ``scipy`` are not installed.
    """
    import os
    rasterio = _require_rasterio()
    from rasterio.crs import CRS

    paths = {
        "dem":   os.path.join(raster_dir, "elevation_cropped.tif"),
        "soil":  os.path.join(raster_dir, "water_storage.tif"),
        "jtemp": os.path.join(raster_dir, "merged_jennings2.tif"),
    }

    for label, p in paths.items():
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"Could not find '{label}' raster at: {p}\n"
                f"Check that raster_dir='{raster_dir}' is correct and that "
                f"the file '{os.path.basename(p)}' exists."
            )

    def _load(path: str, scale: float = 1.0):
        """Read a GeoTIFF; return (array, transform, crs, nodata)."""
        with rasterio.open(path) as src:
            arr = src.read(1).astype(np.float64) * scale
            nodata = src.nodata
            if nodata is not None:
                arr[arr == nodata * scale] = np.nan
            return arr, src.transform, src.crs

    print(f"Loading WBM rasters from '{raster_dir}' ...")

    # ── DEM ─────────────────────────────────────────────────────────────────
    dem_arr, dem_transform, dem_crs = _load(paths["dem"])
    print(f"  DEM loaded        : {dem_arr.shape}  CRS={dem_crs}")

    # ── Soil storage (cm → mm) ───────────────────────────────────────────────
    soil_arr, soil_transform, soil_crs = _load(paths["soil"], scale=10.0)
    print(f"  Soil loaded       : {soil_arr.shape}  CRS={soil_crs}")

    # ── Jennings temperature ─────────────────────────────────────────────────
    jtemp_arr, jtemp_transform, jtemp_crs = _load(paths["jtemp"])
    # The R code sets CRS to 4326 if missing; mirror that here
    if jtemp_crs is None or not jtemp_crs.is_valid:
        from rasterio.crs import CRS as _CRS
        jtemp_crs = _CRS.from_epsg(4326)
        print("  Jennings CRS missing – assumed EPSG:4326")
    print(f"  Jennings loaded   : {jtemp_arr.shape}  CRS={jtemp_crs}")

    # ── Optional reprojection ────────────────────────────────────────────────
    if target_crs is not None:
        dst_crs = CRS.from_user_input(target_crs)

        def _maybe_reproject(arr, transform, src_crs, name, resampling="bilinear"):
            if src_crs != dst_crs:
                print(f"  Reprojecting {name} → {dst_crs.to_string()} ...")
                arr, transform, src_crs = _reproject_raster(
                    arr, transform, src_crs, dst_crs, resampling)
            return arr, transform, src_crs

        dem_arr,   dem_transform,   dem_crs   = _maybe_reproject(dem_arr,   dem_transform,   dem_crs,   "DEM")
        soil_arr,  soil_transform,  soil_crs  = _maybe_reproject(soil_arr,  soil_transform,  soil_crs,  "Soil")
        jtemp_arr, jtemp_transform, jtemp_crs = _maybe_reproject(jtemp_arr, jtemp_transform, jtemp_crs, "Jennings")

    # ── Slope & aspect from DEM (Horn's 8-neighbour method) ──────────────────
    # cell_width / cell_height come from the affine transform.
    # For geographic CRS (degrees) this gives slope in degrees/degree which is
    # dimensionless – same as terra::terrain's behaviour; for projected CRS
    # (metres) both numerator and denominator are in metres.
    cell_w = abs(dem_transform.a)
    cell_h = abs(dem_transform.e)
    print("  Computing slope & aspect (Horn's method) ...")
    slope_arr, aspect_arr = _compute_slope_aspect_horn(dem_arr, cell_w, cell_h)

    # ── Store in module cache ────────────────────────────────────────────────
    _RASTERS.update({
        "dem":            dem_arr,
        "dem_transform":  dem_transform,
        "dem_crs":        dem_crs,
        "slope":          slope_arr,
        "aspect":         aspect_arr,
        "soil":           soil_arr,
        "soil_transform": soil_transform,
        "soil_crs":       soil_crs,
        "jtemp":          jtemp_arr,
        "jtemp_transform":jtemp_transform,
        "jtemp_crs":      jtemp_crs,
    })
    print("  ✓ All rasters loaded and cached.\n")


def _sample_raster(arr: np.ndarray, transform, xs: np.ndarray,
                   ys: np.ndarray) -> np.ndarray:
    """
    Sample raster values at (x, y) point locations using the affine transform.

    Points that fall outside the raster extent return NaN.

    Parameters
    ----------
    arr : np.ndarray  (2-D)
        Raster data array (row 0 = top).
    transform : affine.Affine
        Rasterio affine transform for the raster.
    xs, ys : np.ndarray
        Longitude / latitude coordinates in the same CRS as ``transform``.

    Returns
    -------
    np.ndarray
        Sampled values, one per point.
    """
    from rasterio.transform import rowcol  # lightweight import

    rows, cols = rowcol(transform, xs, ys)
    rows = np.asarray(rows)
    cols = np.asarray(cols)

    height, width = arr.shape
    valid = (rows >= 0) & (rows < height) & (cols >= 0) & (cols < width)

    vals = np.full(len(xs), np.nan)
    vals[valid] = arr[rows[valid], cols[valid]]
    return vals


def extract_point_params(points_df: pd.DataFrame,
                          raster_dir: Optional[str] = None,
                          target_crs: Optional[str] = None) -> pd.DataFrame:
    """
    Extract WBM site parameters from the loaded raster cache.

    If ``load_wbm_rasters()`` has not been called yet, this function calls
    it automatically using the provided (or default) ``raster_dir``.

    Parameters
    ----------
    points_df : pd.DataFrame
        Point locations. Must contain:

        * ``x`` – longitude (decimal degrees, or projected if rasters are
          projected)
        * ``y`` – latitude  (decimal degrees, or projected)

        Any additional columns (e.g. ``GCM``, site IDs) are preserved.

    raster_dir : str or None
        Only used if rasters are not yet loaded. Path to the directory
        holding the three WBM GeoTIFF files. Defaults to
        ``"Data/wbm_rasters"``.

    target_crs : str or None
        Only used if rasters are not yet loaded. Passed to
        ``load_wbm_rasters()``. Leave as ``None`` to keep the native
        EPSG:4326 CRS of the NPS rasters.

    Returns
    -------
    pd.DataFrame
        ``points_df`` with five new columns:

        * ``Elev``    – elevation (m)
        * ``Slope``   – terrain slope (degrees)
        * ``Aspect``  – terrain aspect (degrees, 0–360)
        * ``SWC_Max`` – max soil water storage (mm)
        * ``J_Temp``  – Jennings temperature climatology (°C)

    Notes
    -----
    Slope and aspect are derived from the DEM using Horn's 8-neighbour
    weighted finite-difference algorithm, matching ``terra::terrain``
    with ``neighbors = 8`` in the original R code.

    Soil storage is read from ``water_storage.tif`` and multiplied by 10
    to convert from cm to mm, matching the R code.
    """
    if not _RASTERS:
        rd = raster_dir if raster_dir is not None else _DEFAULT_RASTER_DIR
        load_wbm_rasters(target_crs=target_crs, raster_dir=rd)

    _require_rasterio()   # ensure rasterio is available for rowcol

    xs = points_df["x"].to_numpy(dtype=float)
    ys = points_df["y"].to_numpy(dtype=float)

    elev_vals   = _sample_raster(_RASTERS["dem"],   _RASTERS["dem_transform"],   xs, ys)
    slope_vals  = _sample_raster(_RASTERS["slope"],  _RASTERS["dem_transform"],   xs, ys)
    aspect_vals = _sample_raster(_RASTERS["aspect"], _RASTERS["dem_transform"],   xs, ys)
    soil_vals   = _sample_raster(_RASTERS["soil"],  _RASTERS["soil_transform"],  xs, ys)
    jtemp_vals  = _sample_raster(_RASTERS["jtemp"], _RASTERS["jtemp_transform"], xs, ys)

    # Warn if any points landed outside raster extent
    for name, vals in [("Elev", elev_vals), ("SWC_Max", soil_vals),
                        ("J_Temp", jtemp_vals)]:
        n_nan = np.sum(np.isnan(vals))
        if n_nan > 0:
            warnings.warn(
                f"{n_nan} point(s) returned NaN for '{name}'. "
                "Check that point coordinates fall within the raster extent.",
                stacklevel=2,
            )

    result = points_df.copy()
    result["Elev"]    = elev_vals
    result["Slope"]   = slope_vals
    result["Aspect"]  = aspect_vals
    result["SWC_Max"] = soil_vals
    result["J_Temp"]  = jtemp_vals

     # ── Clamp unrealistic slope values caused by raster artefacts ────────────
    # Slopes ≥ 85° almost always indicate a nodata boundary or water mask edge
    # in the DEM rather than real terrain. Cap at 60° (a very steep but
    # physically plausible hillslope) to prevent the Oudin heat load correction
    # from producing negative PET.
    MAX_REALISTIC_SLOPE = 60.0
    n_bad = (result["Slope"] > MAX_REALISTIC_SLOPE).sum()
    if n_bad > 0:
        warnings.warn(
            f"{n_bad} point(s) had slope > {MAX_REALISTIC_SLOPE}° — likely a "
            f"raster artefact. Clamping to {MAX_REALISTIC_SLOPE}°. "
            f"Affected sites: {result.loc[result['Slope'] > MAX_REALISTIC_SLOPE, 'site'].tolist()}",
            stacklevel=2,
        )
        result["Slope"] = result["Slope"].clip(upper=MAX_REALISTIC_SLOPE)

    return result
