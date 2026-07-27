"""
Precipitation partitioning, snow processes, soil moisture, and storage
reservoir routing for the NPS Water Balance Model.

Ported from the R implementation originally developed by ARC (10/30/2019)
with updates by MGS (4/14/2021) and CSU additions. Based on Dave Thoma's
water balance Excel spreadsheet model.

These functions implement the daily state-update loop (rain/snow split,
snowmelt, snowpack, soil water content, AET, and the CSU linear-storage
addition) that ``wbm.model.nps_wbm`` calls in sequence.
"""

from __future__ import annotations

import math

import numpy as np

# ──────────────────────────────────────────────────────────────────────────────
# 3.  PRECIPITATION PARTITIONING
# ──────────────────────────────────────────────────────────────────────────────

def get_freeze(jtemp: np.ndarray, tmean: np.ndarray) -> np.ndarray:
    """
    Freeze factor (0–1) using Jennings et al. (2018) thresholds.

    0 = all snow, 1 = all rain.

    Parameters
    ----------
    jtemp : array-like
        Jennings temperature climatology (°C).
    tmean : array-like
        Daily mean temperature (°C).
    """
    jtemp = np.asarray(jtemp, dtype=float)
    tmean = np.asarray(tmean, dtype=float)
    lo = jtemp - 3
    hi = jtemp + 3
    freeze = np.where(
        tmean <= lo, 0.0,
        np.where(tmean >= hi, 1.0,
                 (1.0 / (hi - lo)) * (tmean - lo))
    )
    return freeze


def get_rain(ppt: np.ndarray, freeze: np.ndarray) -> np.ndarray:
    """
    Rainfall (mm) = ppt × freeze factor.

    Parameters
    ----------
    ppt : array-like
        Precipitation (mm).
    freeze : array-like
        Freeze factor (0–1).
    """
    return np.asarray(ppt, dtype=float) * np.asarray(freeze, dtype=float)


def get_snow(ppt: np.ndarray, freeze: np.ndarray) -> np.ndarray:
    """
    Snowfall (mm) = ppt × (1 – freeze factor).

    Parameters
    ----------
    ppt : array-like
        Precipitation (mm).
    freeze : array-like
        Freeze factor (0–1).
    """
    return np.asarray(ppt, dtype=float) * (1.0 - np.asarray(freeze, dtype=float))


# ──────────────────────────────────────────────────────────────────────────────
# 4.  SNOW PROCESSES
# ──────────────────────────────────────────────────────────────────────────────

def get_melt(tmean: np.ndarray, jtemp: np.ndarray, hock: float,
             snow: np.ndarray, sp0: float = 0.0) -> np.ndarray:
    """
    Daily snowmelt (mm) using the Hock degree-day approach.

    Melt occurs when tmean ≥ (jtemp − 3) and snowpack > 0.

    Parameters
    ----------
    tmean : array-like
        Daily mean temperature (°C).
    jtemp : array-like
        Jennings temperature climatology (°C); same length as tmean.
    hock : float
        Degree-day melt factor (mm °C⁻¹ day⁻¹).
    snow : array-like
        Daily snowfall (mm).
    sp0 : float, optional
        Initial snowpack (mm). Default 0.

    Returns
    -------
    np.ndarray
        Daily melt (mm).

    References
    ----------
    Hock, R. (2003). J. Hydrology, 282, 104-115.
    https://doi.org/10.1016/S0022-1694(03)00257-9
    """
    tmean = np.asarray(tmean, dtype=float)
    jtemp = np.asarray(jtemp, dtype=float)
    snow = np.asarray(snow, dtype=float)
    n = len(tmean)
    melt = np.zeros(n)
    snowpack = np.zeros(n)

    threshold = jtemp[0] - 3
    potential = (tmean[0] - threshold) * hock
    if tmean[0] < threshold or sp0 == 0:
        melt[0] = 0.0
    else:
        melt[0] = min(potential, sp0)
    snowpack[0] = sp0 + snow[0] - melt[0]

    for i in range(1, n):
        threshold = jtemp[i] - 3
        if tmean[i] < threshold or snowpack[i - 1] == 0:
            melt[i] = 0.0
        else:
            potential = (tmean[i] - threshold) * hock
            melt[i] = min(potential, snowpack[i - 1])
        snowpack[i] = snowpack[i - 1] + snow[i] - melt[i]

    return melt


def get_snowpack(jtemp: np.ndarray, snow: np.ndarray,
                 melt: np.ndarray, sp0: float = 0.0) -> np.ndarray:
    """
    Cumulative snowpack (mm SWE) at each time step.

    Parameters
    ----------
    jtemp : array-like
        Jennings temperature climatology (°C) — kept for API parity with R.
    snow : array-like
        Daily snowfall (mm).
    melt : array-like
        Daily melt (mm).
    sp0 : float, optional
        Initial snowpack (mm). Default 0.

    Returns
    -------
    np.ndarray
        Snowpack at end of each day (mm SWE).
    """
    snow = np.asarray(snow, dtype=float)
    melt = np.asarray(melt, dtype=float)
    n = len(snow)
    snowpack = np.zeros(n)
    sp_i = sp0
    for i in range(n):
        snowpack[i] = sp_i + snow[i] - melt[i]
        sp_i = snowpack[i]
    return snowpack


def get_ablation(tmean: np.ndarray, jtemp: np.ndarray,
                 ab_fac: float, snow: np.ndarray,
                 sp0: float = 0.0) -> np.ndarray:
    """
    Snow ablation (sublimation/vapor loss) from snowpack (mm).

    Water is removed from snowpack but does NOT contribute to liquid runoff.

    Parameters
    ----------
    tmean : array-like
        Daily mean temperature (°C).
    jtemp : array-like
        Jennings temperature climatology (°C).
    ab_fac : float
        Ablation factor (degree-day factor, mm °C⁻¹ day⁻¹).
    snow : array-like
        Daily snowfall (mm).
    sp0 : float, optional
        Initial snowpack (mm). Default 0.

    Returns
    -------
    np.ndarray
        Daily ablation (mm).
    """
    tmean = np.asarray(tmean, dtype=float)
    jtemp = np.asarray(jtemp, dtype=float)
    snow = np.asarray(snow, dtype=float)
    n = len(tmean)
    ablation = np.zeros(n)
    snowpack = np.zeros(n)

    threshold = jtemp[0] - 3
    if tmean[0] < threshold or sp0 == 0:
        ablation[0] = 0.0
    else:
        potential = (tmean[0] - threshold) * ab_fac
        ablation[0] = min(potential, sp0)
    snowpack[0] = sp0 + snow[0] - ablation[0]

    for i in range(1, n):
        threshold = jtemp[i] - 3
        if tmean[i] < threshold or snowpack[i - 1] == 0:
            ablation[i] = 0.0
        else:
            potential = (tmean[i] - threshold) * ab_fac
            ablation[i] = min(potential, snowpack[i - 1])
        snowpack[i] = snowpack[i - 1] + snow[i] - ablation[i]

    return ablation


# ──────────────────────────────────────────────────────────────────────────────
# 6.  SOIL MOISTURE
# ──────────────────────────────────────────────────────────────────────────────

def get_soil(w: np.ndarray, pet: np.ndarray, swc_max: float,
             swc0: float = 0.0) -> np.ndarray:
    """
    Daily soil water content (mm).

    When w > PET the soil is charged; otherwise it drains exponentially.

    Parameters
    ----------
    w : array-like
        Water reaching soil surface (rain + melt) (mm).
    pet : array-like
        Potential evapotranspiration (mm).
    swc_max : float
        Maximum soil water-holding capacity (mm).
    swc0 : float, optional
        Initial soil water content (mm). Default 0.

    Returns
    -------
    np.ndarray
        Soil water content at end of each day (mm).
    """
    w = np.asarray(w, dtype=float)
    pet = np.asarray(pet, dtype=float)
    n = len(w)
    soil = np.zeros(n)
    swc_i = swc0
    w_pet = w - pet

    for i in range(n):
        if w[i] > pet[i]:
            soil[i] = min(w_pet[i] + swc_i, swc_max)
        else:
            soil[i] = swc_i - swc_i * (1 - math.exp(-(pet[i] - w[i]) / swc_max))
        swc_i = soil[i]

    return soil


def get_d_soil(swc: np.ndarray, swc0: float = 0.0) -> np.ndarray:
    """
    Daily change in soil water content (mm).

    Parameters
    ----------
    swc : array-like
        Soil water content time series (mm).
    swc0 : float, optional
        Initial soil water content (mm). Default 0.

    Returns
    -------
    np.ndarray
        Day-to-day change in SWC (mm).
    """
    swc = np.asarray(swc, dtype=float)
    prev = np.concatenate([[swc0], swc[:-1]])
    return swc - prev


def get_aet(w: np.ndarray, pet: np.ndarray,
            swc: np.ndarray, swc0: float = 0.0) -> np.ndarray:
    """
    Actual evapotranspiration (AET) (mm).

    AET = PET when water supply ≥ PET; otherwise AET = w + ΔSWC.

    Parameters
    ----------
    w : array-like
        Water reaching soil surface (rain + melt) (mm).
    pet : array-like
        Potential ET (mm).
    swc : array-like
        Soil water content (mm) — already computed by ``get_soil()``.
    swc0 : float, optional
        Initial soil water content (mm). Default 0.

    Returns
    -------
    np.ndarray
        Daily AET (mm).
    """
    w = np.asarray(w, dtype=float)
    pet = np.asarray(pet, dtype=float)
    swc = np.asarray(swc, dtype=float)
    n = len(w)
    aet = np.zeros(n)
    swc_i = swc0

    for i in range(n):
        if w[i] > pet[i]:
            aet[i] = pet[i]
        else:
            aet[i] = w[i] + swc_i - swc[i]
        swc_i = swc[i]

    return aet


# ──────────────────────────────────────────────────────────────────────────────
# 7.  STORAGE RESERVOIR
# ──────────────────────────────────────────────────────────────────────────────

def get_storage(storage_add: np.ndarray,
                return_rate: float = 1.0) -> dict[str, np.ndarray]:
    """
    Linear storage reservoir: partition water between same-day release and
    carry-over to the next day.

    Parameters
    ----------
    storage_add : array-like
        Daily water added to storage (W − AET − ΔSWC) (mm).
    return_rate : float
        Fraction of total available storage released each day (0–1].
        Default 1.0 (full same-day release → no carry-over).

    Returns
    -------
    dict with keys:
        ``storage_release`` : np.ndarray  (mm, contributes to runoff)
        ``storage_remain``  : np.ndarray  (mm, carry-over to next day)
    """
    storage_add = np.asarray(storage_add, dtype=float)
    n = len(storage_add)
    storage_release = np.zeros(n)
    storage_remain = np.zeros(n)

    storage_release[0] = storage_add[0] * return_rate
    storage_remain[0] = storage_add[0] * (1 - return_rate)

    for i in range(1, n):
        total = storage_remain[i - 1] + storage_add[i]
        storage_release[i] = total * return_rate
        storage_remain[i] = total * (1 - return_rate)

    return {"storage_release": storage_release,
            "storage_remain": storage_remain}
