"""
Quick smoke test for the ``wbm`` package using synthetic data (no rasters
or real climate data required).

Run from the repo root:

    python scripts/smoke_test_wbm.py

This is the same check that used to live as ``_demo()`` at the bottom of
the original monolithic ``nps_wbm_functions.py``, split out here so the
``wbm`` package itself stays import-side-effect-free.
"""

import math

import numpy as np
import pandas as pd

from wbm import nps_wbm


def main():
    rng = np.random.default_rng(42)
    n = 365
    dates = pd.date_range("2000-01-01", periods=n)
    doy = dates.dayofyear.to_numpy(dtype=float)

    tmean = 8 * np.sin(2 * math.pi * (doy - 80) / 365) + 5 + rng.normal(0, 2, n)
    ppt = rng.exponential(3.0, n)

    climate = pd.DataFrame({
        "date":    dates,
        "x":       -105.5,
        "y":       40.0,
        "ppt_mm":  ppt,
        "tmean_C": tmean,
        "GCM":     "historical",
    })

    params = {
        "Elev":    2400,
        "Slope":   10,
        "Aspect":  180,
        "SWC_Max": 150,
        "J_Temp":  1.5,
    }

    result = nps_wbm(climate, params, pet_method="Oudin", to_inches=False)
    print(result[["date", "ppt_mm", "RAIN", "SNOW", "MELT", "PACK",
                  "AET", "RUNOFF", "SOIL"]].tail(10).to_string(index=False))
    print("\nAnnual totals (mm):")
    for col in ["ppt_mm", "RAIN", "SNOW", "MELT", "AET", "RUNOFF"]:
        print(f"  {col}: {result[col].sum():.1f}")


if __name__ == "__main__":
    main()
