"""
Meteorological helpers, daylength, and potential-evapotranspiration (PET) methods.

Ported from the R implementation originally developed by ARC (10/30/2019)
with updates by MGS (4/14/2021) and CSU additions. Based on Dave Thoma's
water balance Excel spreadsheet model.

This module contains the pure meteorology / radiation / PET math used by
``wbm.model.nps_wbm``. None of these functions perform any file I/O, so
they only depend on ``numpy`` and ``pandas``.
"""

from __future__ import annotations

import math
from typing import Optional, Union

import numpy as np
import pandas as pd

# ──────────────────────────────────────────────────────────────────────────────
# 1.  LOW-LEVEL METEOROLOGICAL HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

def get_svp(temp: np.ndarray) -> np.ndarray:
    """
    Saturation Vapor Pressure (kPa).

    Source: eq. 11 in FAO-56.

    Parameters
    ----------
    temp : array-like
        Temperature (°C).

    Returns
    -------
    np.ndarray
        Saturation vapor pressure (kPa).
    """
    temp = np.asarray(temp, dtype=float)
    return 0.6108 * np.exp((17.27 * temp) / (temp + 237.3))


def actual_vp(rhmax: np.ndarray, rhmin: np.ndarray,
              tmax: np.ndarray, tmin: np.ndarray) -> np.ndarray:
    """
    Actual Vapor Pressure (kPa).

    Parameters
    ----------
    rhmax, rhmin : array-like
        Daily max/min relative humidity (%).
    tmax, tmin : array-like
        Daily max/min temperature (°C).

    Returns
    -------
    np.ndarray
        Actual vapor pressure (kPa).
    """
    e_tmax = get_svp(np.asarray(tmax, dtype=float))
    e_tmin = get_svp(np.asarray(tmin, dtype=float))
    return (e_tmin * (rhmax / 100) + e_tmax * (rhmin / 100)) / 2


def atm_press(elev: float) -> float:
    """
    Atmospheric pressure (kPa) estimated from elevation.

    Parameters
    ----------
    elev : float
        Elevation (m).
    """
    return 101.3 * ((293 - 0.0065 * elev) / 293) ** 5.26


def psyc_constant(elev: float) -> float:
    """
    Psychrometric constant (kPa °C⁻¹) from elevation.

    Parameters
    ----------
    elev : float
        Elevation (m).
    """
    return 0.000665 * atm_press(elev)


def vapor_curve(temp: np.ndarray) -> np.ndarray:
    """
    Slope of saturation vapor pressure curve (kPa °C⁻¹).

    Parameters
    ----------
    temp : array-like
        Temperature (°C).
    """
    temp = np.asarray(temp, dtype=float)
    return 4098 * (0.6108 * np.exp(17.27 * temp / (temp + 237.3)) / (temp + 237.3) ** 2)


def clear_sky_rad(doy: np.ndarray, lat: float, elev: float) -> np.ndarray:
    """
    Clear-sky solar radiation (MJ m⁻² day⁻¹).

    Parameters
    ----------
    doy : array-like
        Day-of-year (1–366).
    lat : float
        Latitude (degrees).
    elev : float
        Elevation (m).
    """
    doy = np.asarray(doy, dtype=float)
    d_r = 1 + 0.033 * np.cos((2 * math.pi / 365) * doy)
    declin = 0.409 * np.sin(((2 * math.pi / 365) * doy) - 1.39)
    lat_rad = (math.pi / 180) * lat
    sunset_ang = np.arccos(-np.tan(lat_rad) * np.tan(declin))
    R_a = (
        (24 * 60) / math.pi * 0.0820 * d_r
        * (sunset_ang * np.sin(lat_rad) * np.sin(declin)
           + np.cos(lat_rad) * np.cos(declin) * np.sin(sunset_ang))
    )
    R_so = (0.75 + 2e-5 * elev) * R_a
    return R_so


def outgoing_rad(tmax: np.ndarray, tmin: np.ndarray,
                 R_s: np.ndarray, e_a: np.ndarray,
                 R_so: np.ndarray) -> np.ndarray:
    """
    Net outgoing long-wave radiation (MJ m⁻² day⁻¹).

    Parameters
    ----------
    tmax, tmin : array-like
        Daily max/min temperature (°C).
    R_s : array-like
        Incoming solar radiation (MJ m⁻² day⁻¹).
    e_a : array-like
        Actual vapor pressure (kPa).
    R_so : array-like
        Clear-sky radiation (MJ m⁻² day⁻¹).
    """
    tmax = np.asarray(tmax, dtype=float)
    tmin = np.asarray(tmin, dtype=float)
    R_s = np.asarray(R_s, dtype=float)
    e_a = np.asarray(e_a, dtype=float)
    R_so = np.asarray(R_so, dtype=float)
    R_nl = (
        4.903e-9
        * ((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4) / 2
        * (0.34 - 0.14 * np.sqrt(e_a))
        * (1.35 * (R_s / R_so) - 0.35)
    )
    return R_nl


# ──────────────────────────────────────────────────────────────────────────────
# 2.  DAYLENGTH
# ──────────────────────────────────────────────────────────────────────────────

def get_daylength(dates: pd.Series, lat: Union[float, np.ndarray]) -> np.ndarray:
    """
    Astronomical daylength (hours) from date and latitude.

    Replicates geosphere::daylength() used in the R model.

    Parameters
    ----------
    dates : pd.Series of datetime-like
        Daily dates.
    lat : float or array-like
        Latitude (degrees). Scalar or one value per date.

    Returns
    -------
    np.ndarray
        Daylength in hours for each date.
    """
    dates = pd.to_datetime(dates)
    doy = dates.dt.dayofyear.to_numpy(dtype=float)
    lat = np.asarray(lat, dtype=float)
    if lat.ndim == 0:
        lat = np.full_like(doy, float(lat))

    P = np.arcsin(0.39795 * np.cos(0.2163108 + 2 * np.arctan(0.9671396 * np.tan(0.00860 * (doy - 186)))))
    lat_rad = lat * math.pi / 180
    # clamp argument to [-1, 1] to avoid arccos domain errors
    arg = np.clip((np.sin(0.8333 * math.pi / 180) + np.sin(lat_rad) * np.sin(P))
                  / (np.cos(lat_rad) * np.cos(P)), -1.0, 1.0)
    dayl = 24 - (24 / math.pi) * np.arccos(arg)
    return dayl


# ──────────────────────────────────────────────────────────────────────────────
# 5.  PET METHODS
# ──────────────────────────────────────────────────────────────────────────────

def get_oudin_pet(doy: np.ndarray, lat: float, snowpack: np.ndarray,
                  tmean: np.ndarray, slope: float, aspect: float,
                  shade_coeff: float = 1.0) -> np.ndarray:
    """
    Oudin (2005) daily PET (mm) with topographic heat-load correction.

    PET is zero when snowpack > 2 mm or tmean ≤ −5 °C.

    Parameters
    ----------
    doy : array-like
        Day-of-year (1–366).
    lat : float
        Latitude (degrees).
    snowpack : array-like
        Snowpack (mm SWE).
    tmean : array-like
        Daily mean temperature (°C).
    slope : float
        Site slope (degrees).
    aspect : float
        Site aspect (degrees, 0–360).
    shade_coeff : float, optional
        Shade coefficient (0–1). Default 1.

    Returns
    -------
    np.ndarray
        Daily PET (mm).
    """
    doy = np.asarray(doy, dtype=float)
    snowpack = np.asarray(snowpack, dtype=float)
    tmean = np.asarray(tmean, dtype=float)

    d_r = 1 + 0.033 * np.cos((2 * math.pi / 365) * doy)
    declin = 0.409 * np.sin(((2 * math.pi / 365) * doy) - 1.39)
    lat_rad = (math.pi / 180) * lat
    sunset_ang = np.arccos(np.clip(-np.tan(lat_rad) * np.tan(declin), -1.0, 1.0))
    R_a = (
        (24 * 60) / math.pi * 0.082 * d_r
        * (sunset_ang * np.sin(lat_rad) * np.sin(declin)
           + np.cos(lat_rad) * np.cos(declin) * np.sin(sunset_ang))
    )
    Oudin = np.where(
        snowpack > 2, 0.0,
        np.where(tmean > -5, (R_a * (tmean + 5) * 0.408) / 100, 0.0)
    )

    folded_aspect = abs(180 - abs(aspect - 225))
    lat_rad_scalar = lat * math.pi / 180
    slope_rad = slope * math.pi / 180
    aspect_rad = folded_aspect * math.pi / 180
    heatload = (
        0.339
        + 0.808 * math.cos(lat_rad_scalar) * math.cos(slope_rad)
        - 0.196 * math.sin(lat_rad_scalar) * math.sin(slope_rad)
        - 0.482 * math.cos(aspect_rad) * math.sin(slope_rad)
    )

    #OudinPET = Oudin * heatload * shade_coeff
    OudinPET = np.maximum(Oudin * heatload * shade_coeff, 0.0)
    return OudinPET


def get_hamon_pet(tmean: np.ndarray, daylength: np.ndarray) -> np.ndarray:
    """
    Hamon daily PET (mm).

    Parameters
    ----------
    tmean : array-like
        Daily mean temperature (°C).
    daylength : array-like
        Daylength (hours).

    Returns
    -------
    np.ndarray
        Daily PET (mm).
    """
    tmean = np.asarray(tmean, dtype=float)
    daylength = np.asarray(daylength, dtype=float)
    et_hamon = (
        0.1651
        * (daylength / 12)
        * (216.7 * (6.108 * np.exp((17.26 * tmean) / (tmean + 273.3))))
        / (tmean + 273.3)
    )
    return et_hamon


def get_penman_monteith_pet(tmax: np.ndarray, tmin: np.ndarray,
                             doy: np.ndarray, elev: float, lat: float,
                             rhmax: Optional[np.ndarray] = None,
                             rhmin: Optional[np.ndarray] = None,
                             vp: Optional[np.ndarray] = None,
                             wind: float = 2.2) -> np.ndarray:
    """
    FAO-56 Penman-Monteith daily reference ET (mm).

    Parameters
    ----------
    tmax, tmin : array-like
        Daily max/min temperature (°C).
    doy : array-like
        Day-of-year (1–366).
    elev : float
        Elevation (m).
    lat : float
        Latitude (degrees).
    rhmax, rhmin : array-like, optional
        Daily max/min relative humidity (%). Used for actual VP if ``vp``
        is not supplied.
    vp : array-like, optional
        Daily actual vapor pressure (kPa). Overrides rh-based estimate.
    wind : float, optional
        Mean wind speed (m s⁻¹). Default 2.2 m/s (~5 mph).

    Returns
    -------
    np.ndarray
        Daily reference ET (mm).
    """
    tmax = np.asarray(tmax, dtype=float)
    tmin = np.asarray(tmin, dtype=float)
    tmean = (tmax + tmin) / 2
    doy = np.asarray(doy, dtype=float)

    psyc = psyc_constant(elev)
    delta = vapor_curve(tmean)

    DT = delta / (delta + psyc * (1 + 0.34 * wind))
    PT = psyc / (delta + psyc * (1 + 0.34 * wind))
    TT = (900 / (tmean + 273)) * wind

    e_tmax = get_svp(tmax)
    e_tmin = get_svp(tmin)
    e_s = (e_tmax + e_tmin) / 2

    if vp is not None:
        e_a = np.asarray(vp, dtype=float)
    elif rhmax is not None and rhmin is not None:
        e_a = actual_vp(np.asarray(rhmax), np.asarray(rhmin), tmax, tmin)
    else:
        e_a = e_tmin  # conservative fallback

    R_s = clear_sky_rad(doy, lat, elev)
    R_ns = (1 - 0.23) * R_s
    R_so = clear_sky_rad(doy, lat, elev)
    R_nl = outgoing_rad(tmax, tmin, R_s, e_a, R_so)
    R_n = R_ns - R_nl
    R_ng = 0.408 * R_n

    ET_rad = DT * R_ng
    ET_wind = PT * TT * (e_s - e_a)
    ET_o = ET_rad + ET_wind
    return ET_o


# ──────────────────────────────────────────────────────────────────────────────
# 8.  GROWING DEGREE DAYS  &  DEFICIT
# ──────────────────────────────────────────────────────────────────────────────

def get_gdd(tmean: np.ndarray, t_base: float = 0.0) -> np.ndarray:
    """
    Growing degree days (°C day⁻¹) above a base temperature.

    Parameters
    ----------
    tmean : array-like
        Daily mean temperature (°C).
    t_base : float, optional
        Base temperature (°C). Default 0.

    Returns
    -------
    np.ndarray
        Daily GDD.
    """
    tmean = np.asarray(tmean, dtype=float)
    return np.where(tmean < t_base, 0.0, tmean - t_base)


def get_deficit(pet: np.ndarray, aet: np.ndarray) -> np.ndarray:
    """
    Climatic water deficit = PET − AET (mm).

    Parameters
    ----------
    pet, aet : array-like
        Potential and actual ET (mm).
    """
    return np.asarray(pet, dtype=float) - np.asarray(aet, dtype=float)
