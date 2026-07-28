"""
Incremental per-site caching helpers.

Every expensive per-site step in this repo (GEE downloads in notebook 02,
the Oudin/Penman-Monteith WBM runs in notebooks 03/06, the per-site
PET-multiplier calibration in notebook 04) follows the same shape: loop over
a list of sites, do something slow for each one, concatenate the results
into one CSV.

Before this module existed, resuming that loop after adding new sites meant
an all-or-nothing ``read_prev`` toggle: either skip the step entirely (reuse
the old cache, missing the new sites) or redo every site from scratch
(including ones that were already cached). That's fine for ~35 sites, but
wasteful once the site list grows to ~80+ and only a handful are actually
new.

``load_and_filter_missing`` / ``merge_and_save_cache`` replace that pattern:
figure out which sites are already in the cache, only do the slow work for
the rest, then merge and save once.

Note: this resumes *across* runs (rerunning the notebook after adding new
sites won't redo old ones). It does not checkpoint *within* a run — if the
process is killed partway through the new sites, that batch's progress is
lost and will be redone next time. That's an intentional simplicity
tradeoff; see the docstrings below if you want to add periodic
checkpointing later.
"""

from __future__ import annotations

import os
from typing import Optional, Sequence

import pandas as pd


def load_and_filter_missing(
    all_sites: Sequence[str],
    cache_path: str,
    site_col: str = "site",
    parse_dates: Optional[list] = None,
) -> tuple[Optional[pd.DataFrame], list]:
    """
    Check a per-site cache CSV and return what's missing.

    Parameters
    ----------
    all_sites : sequence of str
        The full list of sites this step needs results for.
    cache_path : str
        Path to the cached CSV (same path you'd pass to ``pd.read_csv`` /
        ``.to_csv`` today).
    site_col : str
        Name of the column identifying each row's site.
    parse_dates : list or None
        Passed through to ``pd.read_csv`` if the cache has date columns.

    Returns
    -------
    (existing_df, missing_sites)
        ``existing_df`` is the full cached DataFrame if the file exists
        (else ``None``). ``missing_sites`` is the subset of ``all_sites``
        not yet present in the cache, in the same order as ``all_sites``.
    """
    if os.path.exists(cache_path):
        existing_df = pd.read_csv(cache_path, parse_dates=parse_dates)
        cached_sites = set(existing_df[site_col].unique())
    else:
        existing_df = None
        cached_sites = set()

    missing = [s for s in all_sites if s not in cached_sites]
    return existing_df, missing


def merge_and_save_cache(
    cache_path: str,
    existing_df: Optional[pd.DataFrame],
    new_rows_df: Optional[pd.DataFrame],
    sort_cols: Optional[list] = None,
) -> pd.DataFrame:
    """
    Combine newly computed rows with whatever was already cached, save, and
    return the combined DataFrame.

    Parameters
    ----------
    cache_path : str
        Where to save the combined result (parent directory is created if
        needed).
    existing_df : DataFrame or None
        Whatever ``load_and_filter_missing`` returned (``None`` if there was
        no prior cache).
    new_rows_df : DataFrame or None
        Freshly computed rows for the sites that were missing. Can be
        ``None``/empty if every site was already cached (nothing new to
        add — this just re-saves ``existing_df`` unchanged).
    sort_cols : list or None
        Columns to sort by before saving (e.g. ``["site", "date"]``).

    Returns
    -------
    DataFrame
        The combined, saved result — assign this back to whatever variable
        name the rest of the notebook expects (e.g. ``climate_gee =
        merge_and_save_cache(...)``).
    """
    frames = [d for d in (existing_df, new_rows_df) if d is not None and len(d) > 0]
    if not frames:
        raise ValueError(
            "Nothing to save: both existing_df and new_rows_df are empty. "
            "If this is the very first run, check that the per-site loop "
            "actually produced rows before calling merge_and_save_cache()."
        )

    combined = pd.concat(frames, ignore_index=True)
    if sort_cols:
        combined = combined.sort_values(sort_cols).reset_index(drop=True)

    out_dir = os.path.dirname(cache_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    combined.to_csv(cache_path, index=False)
    return combined
