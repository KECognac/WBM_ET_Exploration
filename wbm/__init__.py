"""
wbm — NPS Water Balance Model, Python implementation.

Ported from the R implementation originally developed by ARC (10/30/2019)
with updates by MGS (4/14/2021) and CSU additions. Based on Dave Thoma's
water balance Excel spreadsheet model.

This package was split out of the original single-file
``nps_wbm_functions.py`` into task-based modules:

* ``wbm.met_pet``    – meteorological helpers, daylength, and PET methods
                       (Oudin, Hamon, Penman-Monteith), plus GDD/deficit.
* ``wbm.snow_soil``  – precipitation partitioning, snow/melt processes,
                       soil moisture, and the linear storage reservoir.
* ``wbm.model``      – the main daily driver (``nps_wbm``), the
                       multi-point/multi-GCM wrapper (``run_nps_wbm_points``),
                       and the end-to-end convenience pipeline
                       (``run_pipeline``).
* ``wbm.raster_io``  – loading the elevation/soil/Jennings-temperature
                       GeoTIFFs and extracting per-point site parameters
                       (the only module that needs ``rasterio``/``scipy``).

Every public name from those modules is re-exported here, so existing code
that did ``import nps_wbm_functions as wbm`` can be updated to simply
``import wbm`` and continue calling ``wbm.nps_wbm(...)``,
``wbm.load_wbm_rasters(...)``, etc. unchanged.

Typical usage
-------------
    import pandas as pd
    from wbm import nps_wbm, run_nps_wbm_points, load_wbm_rasters, extract_point_params

    load_wbm_rasters(raster_dir="Data/wbm_rasters")
    results = nps_wbm(daily_df=climate_df, point_params=params, pet_method="Oudin")
"""

from .met_pet import (
    get_svp,
    actual_vp,
    atm_press,
    psyc_constant,
    vapor_curve,
    clear_sky_rad,
    outgoing_rad,
    get_daylength,
    get_wind_speed_2m,
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
    get_ablation,
    get_soil,
    get_d_soil,
    get_aet,
    get_storage,
)
from .model import (
    nps_wbm,
    run_nps_wbm_points,
    run_pipeline,
)
from .raster_io import (
    load_wbm_rasters,
    extract_point_params,
)

__all__ = [
    # met_pet
    "get_svp", "actual_vp", "atm_press", "psyc_constant", "vapor_curve",
    "clear_sky_rad", "outgoing_rad", "get_daylength", "get_wind_speed_2m",
    "get_oudin_pet", "get_hamon_pet", "get_penman_monteith_pet",
    "get_gdd", "get_deficit",
    # snow_soil
    "get_freeze", "get_rain", "get_snow", "get_melt", "get_snowpack",
    "get_ablation", "get_soil", "get_d_soil", "get_aet", "get_storage",
    # model
    "nps_wbm", "run_nps_wbm_points", "run_pipeline",
    # raster_io
    "load_wbm_rasters", "extract_point_params",
]

__version__ = "0.1.0"
