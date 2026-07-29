import json, copy

SRC = "10_runoff_implications.ipynb"

with open(SRC) as f:
    nb = json.load(f)

cells = nb["cells"]
by_id = {c["id"]: i for i, c in enumerate(cells)}

def md(text):
    return {"cell_type": "markdown", "id": None, "metadata": {}, "source": text.splitlines(keepends=True)}

def code(text):
    return {"cell_type": "code", "id": None, "metadata": {}, "execution_count": None,
            "outputs": [], "source": text.splitlines(keepends=True)}

import uuid
def new_id():
    return uuid.uuid4().hex[:8]

def set_cell(cell_id, new_cell):
    idx = by_id[cell_id]
    new_cell["id"] = cell_id
    new_cell["metadata"] = cells[idx].get("metadata", {})
    cells[idx] = new_cell

# ---------------------------------------------------------------------------
# Cell df5495d9 -- intro markdown: mention six variants + new seasonal section
# ---------------------------------------------------------------------------
intro = '''## Notebook 10 — Water-balance implications of AET method choice

**Why this matters:** in the NPS WBM's bucket structure, actual evapotranspiration
(AET) is the model's primary *outflow* competing with runoff generation. Concretely
(see `wbm/snow_soil.py`): each day, `storage_add = W - AET - dSOIL` (where `W` is rain
+ snowmelt reaching the soil), and `RUNOFF = storage_release(storage_add) + direct_frac
* RAIN`. Higher AET means less water left over for storage release, so **at the
annual-to-multiyear scale, underestimating AET directly inflates the WBM's runoff
estimate** (soil-moisture and snowpack carry-over roughly net to zero over several
years, so `RUNOFF ≈ PPT - AET` is a good long-run approximation, and this notebook
checks that against the model's own literal `RUNOFF` output rather than assuming it).

Notebooks 04/07/08 established that **default (uncalibrated) Oudin systematically
underestimates AET** relative to OpenET and flux towers, and that per-site- and
per-ecosystem-calibrated Oudin and Penman-Monteith (default and Kc-calibrated) all
correct much of that bias to varying degrees. This notebook asks the practical
follow-on question, designed around the same **six-variant** comparison used in
notebook 08: **if the WBM is currently being run with uncalibrated Oudin, how much is
its runoff estimate being inflated as a result — and how much of that gets corrected
by each alternative?**

**A note on the two Penman-Monteith Kc-calibrated variants:** notebook 07's per-site
calibration cache is currently stale — its incremental caching only checks whether a
site's *name* is already present, not whether the climate data it was fit against has
changed, so it silently skipped re-fitting `Kc` after the gridMET wind-data refetch
earlier in this project. The load cell below detects this automatically (by comparing
file timestamps) and excludes `pmcal`/`pmecocal` from every comparison until notebook
07 is re-run with its stale caches deleted first — at that point, just re-running this
notebook will pick the fresh data up with no further edits needed.

A second question, addressed in the new "Seasonality" section below: **is that
runoff-overestimation bias spread evenly across the year, or concentrated in
particular months?** AET's share of the water balance varies seasonally (e.g. low in
winter when little water is available for ET, high in the growing season), so a
temperature-based method's AET bias plausibly does too — and if the biggest bias
lines up with the months that already generate the most runoff (e.g. snowmelt), the
annual-total framing above could be masking a much larger practical error at the
specific times of year that matter most for downstream water availability.

This is a read-only analysis: it reuses the full 2016–2023 daily WBM outputs already
cached by notebooks 03 (default Oudin), 04 (per-site/per-ecosystem-calibrated Oudin),
06 (Penman-Monteith), and 07 (per-site/per-ecosystem Kc-calibrated Penman-Monteith) —
no new WBM runs, no GEE access. `RUNOFF` is already a per-day output column in each of
those caches, so this notebook just aggregates and compares it directly, rather than
re-deriving it from `PPT - AET`.

**Prerequisites:** run notebooks 01→02→03 at least once (default Oudin WBM cache),
notebook 04 at least once (`wbm_results_site_cal_{OBJECTIVE}.csv` and
`wbm_results_ecosystem_cal_{OBJECTIVE}.csv`), notebook 06 at least once
(`wbm_results_penman_monteith.csv`), and notebook 07 at least once
(`wbm_results_site_cal_pm_{OBJECTIVE}.csv` and
`wbm_results_ecosystem_cal_pm_{OBJECTIVE}.csv`).'''
set_cell("df5495d9", md(intro))

# ---------------------------------------------------------------------------
# Cell ea990ddf -- load all six variants
# ---------------------------------------------------------------------------
load_code = '''IN_TO_MM = 25.4
COLS = ["date", "site", "ecosystem", "state", "ppt_mm", "AET", "RUNOFF"]

with open("../Data/gridmet_cache/last_objective.txt") as f:
    CAL_OBJECTIVE = f.read().strip()
with open("../Data/gridmet_cache/last_objective_pm.txt") as f:
    CAL_OBJECTIVE_PM = f.read().strip()
if CAL_OBJECTIVE_PM != CAL_OBJECTIVE:
    print(f"WARNING -- Oudin calibration objective ('{CAL_OBJECTIVE}') and Penman-Monteith "
          f"Kc calibration objective ('{CAL_OBJECTIVE_PM}') differ -- the calibrated "
          "variants below are not using the same objective function (same check as "
          "notebook 08's banner).")

# Six variants, matching notebook 08's comparison exactly.
VARIANT_FILES = {
    "default":  "wbm_results_flux_towers_2016_2023.csv",
    "sitecal":  f"wbm_results_site_cal_{CAL_OBJECTIVE}.csv",
    "ecocal":   f"wbm_results_ecosystem_cal_{CAL_OBJECTIVE}.csv",
    "pm":       "wbm_results_penman_monteith.csv",
    "pmcal":    f"wbm_results_site_cal_pm_{CAL_OBJECTIVE_PM}.csv",
    "pmecocal": f"wbm_results_ecosystem_cal_pm_{CAL_OBJECTIVE_PM}.csv",
}
VARIANT_LABELS = {
    "default":  "Oudin (default)",
    "sitecal":  f"Oudin (per-site calibrated, {CAL_OBJECTIVE})",
    "ecocal":   f"Oudin (per-ecosystem calibrated, {CAL_OBJECTIVE})",
    "pm":       "Penman-Monteith (default)",
    "pmcal":    f"Penman-Monteith (per-site Kc calibrated, {CAL_OBJECTIVE_PM})",
    "pmecocal": f"Penman-Monteith (per-ecosystem Kc calibrated, {CAL_OBJECTIVE_PM})",
}
ALT_COLORS = {
    "sitecal":  "seagreen",
    "ecocal":   "goldenrod",
    "pm":       "tomato",
    "pmcal":    "mediumpurple",
    "pmecocal": "slategrey",
}

# All caches store fluxes in inches (the repo's default TO_INCHES = True), despite the
# _mm-suffixed column names -- converted to mm here (x 25.4) to match units used
# elsewhere in this project's reports.
#
# Auto-detect missing/stale PM-Kc-calibrated caches: Oudin never uses wind, so its
# calibrated variants (sitecal/ecocal) are safe regardless of when they were last
# (re)run. Penman-Monteith's Kc calibration (notebook 07) DOES depend on wind -- and
# its incremental caching only checks site PRESENCE, not whether the underlying
# climate data changed, so it silently skipped re-fitting after the gridMET wind-data
# refetch earlier in this project. Two possible states here: (a) the stale cache
# still exists (file timestamp predates the climate refetch) or (b) it's been deleted
# in preparation for re-running notebook 07 (file doesn't exist yet). Either way,
# exclude it from every comparison below until a fresh cache appears.
CLIMATE_PATH = "../Data/gridmet_cache/climate_gee_flux_towers_2016_2023.csv"
_climate_mtime = os.path.getmtime(CLIMATE_PATH) if os.path.exists(CLIMATE_PATH) else 0

wbm = {}
ALTERNATIVES = ["sitecal", "ecocal", "pm"]  # always trustworthy (Oudin has no wind dependency)
for key, fname in VARIANT_FILES.items():
    fpath = f"../Data/gridmet_cache/{fname}"

    if key in ("pmcal", "pmecocal"):
        if not os.path.exists(fpath):
            print(f"\\nNOTE -- '{VARIANT_LABELS[key]}' cache ({fname}) not found -- looks "
                  "like it's been deleted in preparation for re-running notebook 07. Excluded "
                  "from every comparison below until a fresh cache is generated; just re-run "
                  "this notebook afterward and it will be picked up automatically.")
            continue
        if os.path.getmtime(fpath) < _climate_mtime:
            print(f"\\nWARNING -- '{VARIANT_LABELS[key]}' cache ({fname}) predates the "
                  f"current gridMET climate data ({CLIMATE_PATH}) -- notebook 07's "
                  "incremental per-site caching silently skipped re-fitting Kc after the "
                  "wind-data refetch (it only checks site presence, not whether the climate "
                  "data changed). EXCLUDED from every comparison below until notebook 07 is "
                  "re-run with its stale caches (kc_site_*.csv, kc_ecosystem_*.csv, "
                  "wbm_results_(site|ecosystem)_cal_pm_*.csv, "
                  "wbm_monthly_(site|ecosystem)_cal_pm_*.csv) deleted first.")
            continue
        ALTERNATIVES.append(key)

    df = pd.read_csv(fpath, usecols=COLS, parse_dates=["date"])
    for c in ("ppt_mm", "AET", "RUNOFF"):
        df[c] = df[c] * IN_TO_MM
    wbm[key] = df
    print(f"{key:10s} ({VARIANT_LABELS[key]:58s}): {df.shape}")

print(f"\\nVariants included in this run's comparisons: {ALTERNATIVES}")'''
set_cell("ea990ddf", code(load_code))

# ---------------------------------------------------------------------------
# Cell b33e1498 -- mean annual balance, all six variants
# ---------------------------------------------------------------------------
annual_code = '''def mean_annual_balance(df, label):
    """Calendar-year sums per site, then averaged across years -> one row per site."""
    yearly = (
        df.assign(year=df["date"].dt.year)
        .groupby(["site", "year"], as_index=False)[["ppt_mm", "AET", "RUNOFF"]]
        .sum()
    )
    out = (
        yearly.groupby("site", as_index=False)[["ppt_mm", "AET", "RUNOFF"]]
        .mean()
        .rename(columns={
            "ppt_mm": f"ppt_{label}", "AET": f"aet_{label}", "RUNOFF": f"runoff_{label}",
        })
    )
    return out


ann = {key: mean_annual_balance(df, key) for key, df in wbm.items()}
eco_lookup = wbm["default"][["site", "ecosystem", "state"]].drop_duplicates()

balance = ann["default"]
for key in ALTERNATIVES:
    balance = balance.merge(ann[key], on="site", how="inner")
balance = balance.merge(eco_lookup, on="site", how="left")

# Sanity check: RUNOFF ~ PPT - AET at the multi-year mean scale (storage/snowpack
# carry-over should mostly net out over 8 years) -- this is NOT assumed anywhere
# above, just checked here against the model's own literal RUNOFF output. Only
# checked for variants actually included below (ALTERNATIVES already excludes any
# stale PM-Kc cache -- see the WARNING printed in the load cell above).
for label in ["default"] + ALTERNATIVES:
    resid = (balance[f"ppt_{label}"] - balance[f"aet_{label}"]) - balance[f"runoff_{label}"]
    print(f"{label:10s}: PPT - AET vs RUNOFF, mean residual = {resid.mean():+.1f} mm/yr "
          f"(median abs residual = {resid.abs().median():.1f} mm/yr, n={len(balance)} sites)")

print(f"\\nSites with all {1 + len(ALTERNATIVES)} included variants matched: {len(balance)}")
balance.to_csv("../Data/runoff_implications/annual_water_balance_by_site.csv", index=False)'''
set_cell("b33e1498", code(annual_code))

# ---------------------------------------------------------------------------
# Cell 830f7201 -- delta_runoff / overestimate %, all five alternatives
# ---------------------------------------------------------------------------
delta_code = '''# "% of the alternative's own runoff estimate" is a natural "how far off is the status
# quo" framing, but blows up at sites where the alternative's runoff is itself near
# zero (e.g. very dry sites) -- restrict that specific ratio to sites with a
# meaningfully nonzero denominator (>= 10 mm/yr) so a handful of near-zero-runoff sites
# don't dominate the mean. The %-of-precipitation version below has no such issue
# (precipitation is never near zero) and is reported for all sites.
RUNOFF_FLOOR = 10.0  # mm/yr

for alt in ALTERNATIVES:
    balance[f"delta_runoff_{alt}"] = balance["runoff_default"] - balance[f"runoff_{alt}"]
    balance[f"pct_overest_{alt}_of_runoff"] = np.where(
        balance[f"runoff_{alt}"] >= RUNOFF_FLOOR,
        100 * balance[f"delta_runoff_{alt}"] / balance[f"runoff_{alt}"],
        np.nan,
    )
    balance[f"pct_overest_{alt}_of_ppt"] = (
        100 * balance[f"delta_runoff_{alt}"] / balance["ppt_default"].replace(0, np.nan)
    )

print("── Nationwide: how much does uncalibrated-Oudin runoff shrink if switched? ──")
for alt in ALTERNATIVES:
    label = VARIANT_LABELS[alt]
    d_mm   = balance[f"delta_runoff_{alt}"]
    d_pctR = balance[f"pct_overest_{alt}_of_runoff"].dropna()
    d_pctP = balance[f"pct_overest_{alt}_of_ppt"]
    n_excluded = balance[f"runoff_{alt}"].lt(RUNOFF_FLOOR).sum()
    print(f"\\nSwitch default Oudin -> {label}:")
    print(f"  mean  Δrunoff = {d_mm.mean():+7.1f} mm/yr   "
          f"median Δrunoff = {d_mm.median():+7.1f} mm/yr")
    print(f"  mean  overestimate = {d_pctR.mean():+6.1f}% of the {label} runoff estimate  "
          f"(median = {d_pctR.median():+6.1f}%, {n_excluded} near-zero-runoff site(s) excluded)")
    print(f"  mean  overestimate = {d_pctP.mean():+6.1f}% of mean annual precip  "
          f"(median = {d_pctP.median():+6.1f}%)")

balance.to_csv("../Data/runoff_implications/annual_water_balance_by_site.csv", index=False)
print("\\nSaved to Data/runoff_implications/annual_water_balance_by_site.csv")'''
set_cell("830f7201", code(delta_code))

# ---------------------------------------------------------------------------
# Cell 964c8290 -- by-ecosystem bar chart, five bars
# ---------------------------------------------------------------------------
eco_bar_code = '''_included = ["default"] + ALTERNATIVES  # only variants actually merged into `balance`
eco_cols = (
    ["ppt_default"]
    + [f"aet_{k}" for k in _included]
    + [f"runoff_{k}" for k in _included]
    + [f"delta_runoff_{a}" for a in ALTERNATIVES]
    + [f"pct_overest_{a}_of_ppt" for a in ALTERNATIVES]
)
eco_summary = (
    balance.groupby("ecosystem")[eco_cols].mean().round(1)
    .join(balance.groupby("ecosystem").size().rename("n_sites"))
    .reset_index()
    .sort_values("n_sites", ascending=False)
)
print(eco_summary.to_string(index=False))
eco_summary.to_csv("../Data/runoff_implications/annual_water_balance_by_ecosystem.csv", index=False)

x = np.arange(len(eco_summary))
w = 0.15
fig, ax = plt.subplots(figsize=(10.5, 4.5), constrained_layout=True)
offsets = np.linspace(-2, 2, len(ALTERNATIVES)) * w
for alt, off in zip(ALTERNATIVES, offsets):
    ax.bar(x + off, eco_summary[f"pct_overest_{alt}_of_ppt"], width=w,
           color=ALT_COLORS[alt], label=f"Switch to {VARIANT_LABELS[alt]}")
ax.axhline(0, color="grey", linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels(eco_summary["ecosystem"], rotation=30, ha="right", fontsize=8.5)
ax.set_ylabel("Mean runoff overestimate from using\\nuncalibrated Oudin (% of annual precip)", fontsize=9)
ax.set_title("How much does uncalibrated Oudin overestimate annual runoff?\\n"
             "by ecosystem, as % of mean annual precipitation (six-variant comparison)",
             fontsize=10.5, fontweight="bold")
ax.legend(fontsize=7, frameon=False, ncol=2)
ax.spines[["top", "right"]].set_visible(False)
plt.savefig("../Data/runoff_implications/runoff_overestimate_by_ecosystem.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved to Data/runoff_implications/runoff_overestimate_by_ecosystem.png")'''
set_cell("964c8290", code(eco_bar_code))

# ---------------------------------------------------------------------------
# Cell 193095ea -- site-level scatter, grid sized to however many
# alternatives are actually included this run (3-5 depending on whether the
# PM-Kc variants passed the staleness check above)
# ---------------------------------------------------------------------------
scatter_code = '''ecosystems_all = sorted(balance["ecosystem"].dropna().unique())
tab_colors = plt.cm.tab10.colors
eco_colors = {e: tab_colors[i % 10] for i, e in enumerate(ecosystems_all)}

n_alt = len(ALTERNATIVES)
ncols = min(n_alt, 3)
nrows = -(-n_alt // ncols)  # ceil division
fig, axes = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 4.2 * nrows), constrained_layout=True)
axes = np.atleast_1d(axes).ravel()
panels = [(f"runoff_{alt}", VARIANT_LABELS[alt]) for alt in ALTERNATIVES]

for ax, (col, label) in zip(axes, panels):
    for eco, grp in balance.groupby("ecosystem"):
        ax.scatter(grp[col], grp["runoff_default"], s=40, color=eco_colors.get(eco, "grey"),
                   edgecolors="k", linewidths=0.4, alpha=0.85, label=eco, zorder=3)
    lims = [0, max(balance[col].max(), balance["runoff_default"].max()) * 1.08]
    ax.plot(lims, lims, color="grey", linestyle="--", linewidth=1, zorder=1)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(f"{label} annual runoff (mm/yr)", fontsize=8.5)
    ax.set_ylabel("Uncalibrated-Oudin annual runoff (mm/yr)", fontsize=8.5)
    ax.set_title(f"vs. {label}", fontsize=9.5, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
for ax in axes[n_alt:]:
    ax.set_visible(False)

handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=7,
                      markerfacecolor=eco_colors[e], markeredgecolor="k", label=e)
           for e in ecosystems_all]
fig.legend(handles=handles, loc="lower right", ncol=1,
           bbox_to_anchor=(0.98, 0.02), fontsize=8, frameon=False)
fig.suptitle("Points above the 1:1 line = uncalibrated Oudin overestimates annual runoff at that site",
             fontsize=11.5, fontweight="bold")
plt.savefig("../Data/runoff_implications/runoff_scatter_default_vs_alternatives.png",
            dpi=150, bbox_inches="tight")
plt.show()
print("Saved to Data/runoff_implications/runoff_scatter_default_vs_alternatives.png")

for alt in ALTERNATIVES:
    n_above = (balance["runoff_default"] > balance[f"runoff_{alt}"]).sum()
    print(f"Sites where uncalibrated Oudin overestimates runoff vs. {VARIANT_LABELS[alt]}: "
          f"{n_above} / {len(balance)}")'''
set_cell("193095ea", code(scatter_code))

# ---------------------------------------------------------------------------
# NEW cells -- seasonality section, inserted before the Takeaways cell
# ---------------------------------------------------------------------------
seasonality_intro_md = '''### Seasonality of the runoff-overestimation bias

The annual-total comparison above shows *how much* switching methods would change
mean annual runoff, but says nothing about *when in the year* that correction would
actually apply. If uncalibrated Oudin's AET underestimation is concentrated in the
months that already produce the most runoff (e.g. spring snowmelt), the practical
consequence is larger than the annual-average framing suggests — a flood-frequency or
peak-flow estimate driven by those specific months would inherit most of the bias,
even though the *annual total* correction looks the same either way. Conversely, if
the bias is spread evenly across the year, the annual-total number already tells the
whole story.

This section repeats the same `RUNOFF`/`AET` comparison at **calendar-month**
resolution: for each site, average each calendar month's total `PPT`/`AET`/`RUNOFF`
across all 8 years (a monthly climatology), then look at how the runoff-overestimation
bias (`Δrunoff = runoff_default − runoff_alt`) and the AET bias
(`ΔAET = AET_alt − AET_default`) vary month to month, nationwide and by
ecosystem.'''

monthly_climatology_code = '''def monthly_climatology(df, label):
    """Calendar-month totals per site per year, then averaged across years ->
    one row per site per month (1-12) -- the 'typical monthly water balance'
    each variant implies at each site."""
    monthly_totals = (
        df.assign(year=df["date"].dt.year, month=df["date"].dt.month)
        .groupby(["site", "year", "month"], as_index=False)[["ppt_mm", "AET", "RUNOFF"]]
        .sum()
    )
    out = (
        monthly_totals.groupby(["site", "month"], as_index=False)[["ppt_mm", "AET", "RUNOFF"]]
        .mean()
        .rename(columns={
            "ppt_mm": f"ppt_{label}", "AET": f"aet_{label}", "RUNOFF": f"runoff_{label}",
        })
    )
    return out


mon = {key: monthly_climatology(df, key) for key, df in wbm.items()}

monthly = mon["default"]
for key in ALTERNATIVES:
    monthly = monthly.merge(mon[key], on=["site", "month"], how="inner")
monthly = monthly.merge(eco_lookup, on="site", how="left")

for alt in ALTERNATIVES:
    monthly[f"delta_runoff_{alt}"] = monthly["runoff_default"] - monthly[f"runoff_{alt}"]
    monthly[f"delta_aet_{alt}"]    = monthly[f"aet_{alt}"] - monthly["aet_default"]

print(f"Monthly climatology rows: {monthly.shape}  "
      f"({monthly['site'].nunique()} sites x 12 months, should be sites x 12)")
monthly.to_csv("../Data/runoff_implications/monthly_water_balance_by_site.csv", index=False)
print("Saved to Data/runoff_implications/monthly_water_balance_by_site.csv")'''

monthly_nationwide_code = '''MONTH_NAMES = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

nation_monthly = monthly.groupby("month")[
    ["runoff_default", "aet_default"]
    + [f"runoff_{a}" for a in ALTERNATIVES]
    + [f"aet_{a}" for a in ALTERNATIVES]
    + [f"delta_runoff_{a}" for a in ALTERNATIVES]
    + [f"delta_aet_{a}" for a in ALTERNATIVES]
].mean()
nation_monthly.to_csv("../Data/runoff_implications/monthly_water_balance_nationwide.csv")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)

# Panel 1: monthly runoff-overestimation bias, all five alternatives.
ax = axes[0]
for alt in ALTERNATIVES:
    ax.plot(nation_monthly.index, nation_monthly[f"delta_runoff_{alt}"],
             marker="o", markersize=4, color=ALT_COLORS[alt], label=VARIANT_LABELS[alt])
ax.axhline(0, color="grey", linewidth=0.8)
ax.set_xticks(range(1, 13)); ax.set_xticklabels(MONTH_NAMES, fontsize=8)
ax.set_ylabel("Mean Δrunoff = runoff_default - runoff_alt (mm/month)", fontsize=9)
ax.set_title("Monthly runoff overestimate from uncalibrated Oudin\\n(nationwide mean across sites)",
             fontsize=10, fontweight="bold")
ax.legend(fontsize=6.3, frameon=False)
ax.spines[["top", "right"]].set_visible(False)

# Panel 2: bias overlaid on the alternative's own monthly runoff magnitude, to
# visually test whether the bias tracks the size of runoff itself.
ax = axes[1]
ax2 = ax.twinx()
ax2.bar(nation_monthly.index, nation_monthly["runoff_default"], color="lightsteelblue",
        alpha=0.6, width=0.7, zorder=1, label="Uncalibrated-Oudin monthly runoff (right axis)")
for alt in ("sitecal", "pm"):
    ax.plot(nation_monthly.index, nation_monthly[f"delta_runoff_{alt}"],
             marker="o", markersize=4, color=ALT_COLORS[alt], label=VARIANT_LABELS[alt], zorder=3)
ax.axhline(0, color="grey", linewidth=0.8, zorder=2)
ax.set_xticks(range(1, 13)); ax.set_xticklabels(MONTH_NAMES, fontsize=8)
ax.set_ylabel("Mean Δrunoff (mm/month)", fontsize=9)
ax2.set_ylabel("Mean monthly runoff, uncal. Oudin (mm)", fontsize=8, color="steelblue")
ax.set_title("Runoff-overestimation bias vs. runoff magnitude, by month", fontsize=10, fontweight="bold")
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, fontsize=6.3, frameon=False, loc="upper left")
ax.spines[["top"]].set_visible(False)

plt.savefig("../Data/runoff_implications/monthly_runoff_bias_nationwide.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved to Data/runoff_implications/monthly_runoff_bias_nationwide.png")

# ── Directly test the motivating hypothesis: is the overestimation bias
# biggest in the months that already generate the most runoff? ─────────────
print("\\nAcross the 12 calendar months (nationwide mean), correlation between "
      "uncalibrated-Oudin's\\nown monthly runoff magnitude and the size of its "
      "overestimation bias vs. each alternative:")
for alt in ALTERNATIVES:
    corr = nation_monthly["runoff_default"].corr(nation_monthly[f"delta_runoff_{alt}"])
    print(f"  runoff_default vs delta_runoff_{alt:10s}: r = {corr:+.2f}")

peak_runoff_month = MONTH_NAMES[int(nation_monthly["runoff_default"].idxmax()) - 1]
print(f"\\nMonth with largest mean runoff (uncal. Oudin): {peak_runoff_month} "
      f"({nation_monthly['runoff_default'].max():.1f} mm/month)")
for alt in ALTERNATIVES:
    worst_m = nation_monthly[f"delta_runoff_{alt}"].idxmax()
    print(f"Month with largest Δrunoff vs {VARIANT_LABELS[alt]}: "
          f"{MONTH_NAMES[int(worst_m) - 1]} ({nation_monthly[f'delta_runoff_{alt}'].max():+.1f} mm/month)")

print("\\nFor context, monthly AET bias (ΔAET = AET_alt - AET_default, positive means "
      "the alternative\\nadds back AET the default was missing) across the same months:")
for alt in ALTERNATIVES:
    worst_m = nation_monthly[f"delta_aet_{alt}"].idxmax()
    print(f"  Largest ΔAET vs {VARIANT_LABELS[alt]:58s}: "
          f"{MONTH_NAMES[int(worst_m) - 1]} ({nation_monthly[f'delta_aet_{alt}'].max():+.1f} mm/month)")'''

monthly_heatmap_code = '''eco_month = (
    monthly.groupby(["ecosystem", "month"])[[f"delta_runoff_{a}" for a in ALTERNATIVES]]
    .mean()
    .reset_index()
)
eco_month.to_csv("../Data/runoff_implications/monthly_water_balance_by_ecosystem.csv", index=False)

eco_order = eco_summary["ecosystem"].tolist()  # sorted by n_sites desc, from the annual section

fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
for ax, alt in zip(axes, ("sitecal", "pm")):
    pivot = (
        eco_month.pivot(index="ecosystem", columns="month", values=f"delta_runoff_{alt}")
        .reindex(eco_order)
    )
    vmax = max(abs(np.nanmin(pivot.values)), abs(np.nanmax(pivot.values)), 1)
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(12)); ax.set_xticklabels(MONTH_NAMES, fontsize=8)
    ax.set_yticks(range(len(eco_order))); ax.set_yticklabels(eco_order, fontsize=8)
    ax.set_title(f"Δrunoff (default − {VARIANT_LABELS[alt]})\\nby ecosystem x month (mm)",
                 fontsize=9.5, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label("mm/month", fontsize=7.5)
plt.savefig("../Data/runoff_implications/monthly_bias_heatmap_by_ecosystem.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved to Data/runoff_implications/monthly_bias_heatmap_by_ecosystem.png")
print("\\nRed = uncalibrated Oudin overestimates runoff that month (positive Δrunoff); "
      "blue = it underestimates it (negative Δrunoff) relative to that alternative.")'''

new_cells = [
    md(seasonality_intro_md),
    code(monthly_climatology_code),
    code(monthly_nationwide_code),
    code(monthly_heatmap_code),
]
for c in new_cells:
    c["id"] = new_id()

takeaways_idx = by_id["a0804a19"]
for i, c in enumerate(new_cells):
    cells.insert(takeaways_idx + i, c)

nb["cells"] = cells

with open(SRC, "w") as f:
    json.dump(nb, f, indent=1)
    f.write("\n")

print("Wrote", SRC, "with", len(cells), "cells")
