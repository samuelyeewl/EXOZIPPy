# Model files: EXOFASTv2-compatible model and residual text files (plan)

Status (2026-10-14): IMPLEMENTED. The decisions are in section 8; sections 4, 6 and 9
record where the build differs from the plan.

**Summary.** EXOFASTv2 writes the model, residuals and "pretty" model of every dataset as
plain text, so users can make their own figures without rebuilding the model. EXOZIPPy
writes only PDFs, the trace and `<prefix>_results.csv`. This plan adds the same files to
EXOZIPPy, **with EXOFASTv2's directory, file names, columns, units, number formats and
index padding**, so a script written against EXOFASTv2 output reads EXOZIPPy output
unchanged.

They are written at the build start, at the polished start and at the median draw of the
reported posterior. A console script also writes them from a finished fit's saved trace.

---

## 1. The reference: which EXOFASTv2

EXOFASTv2's model files differ between branches, so the reference has to be pinned:

| | `master` (2022-07) | `samsmods` fork (f419edcd) | `linkpars_current` (60909ec9, 2026-01) |
|---|---|---|---|
| directory | beside the PDFs | beside the PDFs | **`modelfiles/` subdirectory** |
| dataset index | `transit_%03d`, `telescope_%02d` | unpadded (`transit_0`) | `transit_%03d`, `telescope_%02d` |
| `detrendedmodel.telescope_*`, `trend.telescope_*` | no | no | yes |
| `residuals.telescope_*` error column | raw error | `sqrt(err^2 + jitter variance)` | raw error |
| `sed.residuals.txt` columns | -- | 7 (`mistsed.pro`) | 8 (`mistmultised.pro`, section 3.4) |
| `atmosphere.%03d.txt` | no | no | yes |
| `mist.track` | no | yes (3 columns) | no (`massradius_mist` writes 4 columns, but `exofast_chi2v2` never asks it to) |

**The reference is `linkpars_current`**, the current upstream line and the one the
dilution work followed (`docs/transit_dilution_plan.md`). It is also the only branch that
writes to `modelfiles/`. Two `samsmods` outputs have no upstream equivalent and are kept
because the sample relies on them. The first is `mist.track`, in the 4-column format
upstream's own writer uses (section 3.5). The second is the jitter-inflated RV error, as
an extra fourth column after upstream's raw error (section 3.2). Line numbers below are
`linkpars_current`'s.

## 2. Layout, names, formats

**Directory.** `<dirname(prefix)>/modelfiles/`, created when missing, as
`exofastv2.pro:1337` does. A prefix with no directory part gives `./modelfiles/`.

**Names.** EXOFASTv2's base is `modelfiles/<basename(prefix)><tag>`, where its prefix ends
in `.`. EXOZIPPy's prefix does not end in a dot, so the base is
`modelfiles/<basename(prefix)>.<tag>`, and every file is `<base>.<kind>.<dataset>.txt`, for
example `fitresults/modelfiles/toi2876.mcmc.model.transit_000.txt`.

| tag | point | EXOFASTv2's tag |
|---|---|---|
| `start` | `build_internal_start` (the seeds the user wrote) | `start` |
| `polished` | `internal_start` (the seed polish's result), only when the polish ran | `amoeba` (named `polished` here, decision 1) |
| `mcmc` | the median draw of the reported posterior (section 5) | `mcmc` (its best-likelihood draw) |
| `mcmc_<mode_suffix(k)>` | the median draw of mode k, when the posterior is multimodal | -- |

**Dataset numbers** are the component's own element order: `transit_000` is
`transit.names[0]`, `telescope_00` is the first RV instrument, and `planet_00` is
`system.planet` element 0. EXOFASTv2 numbers files in its own (glob) order, so the numbers
need not match an EXOFASTv2 fit of the same data. Those files have no header, so EXOZIPPy
also writes `modelfiles/<basename(prefix)>.index.txt`, an addition, not an EXOFASTv2
file. It has one line per numbered dataset: `transit_000  <name>  <data file path>`.
The run log carries the same mapping, as EXOFASTv2's does.

**Formats.** These are IDL format codes translated literally: `f0.N` -> `%.Nf`, `e0.N` ->
`%.Ne`, columns separated by one space (`x`). Files EXOFASTv2 writes with `/nocomment` get
**no header line**. The two files that have a comment line get exactly EXOFASTv2's.

**Units are EXOFASTv2's, whatever the input files used.**
- Times are BJD_TDB: EXOZIPPy's internal times, after any `time_offset:` /
  `time_scale:` / `time_frame:` conversion. EXOFASTv2 requires BJD_TDB input, so this is
  the same quantity.
- Transit flux is in the data file's own (normalized) flux units, as EXOFASTv2 writes
  `transit.flux`.
- RVs are in m/s, converted from the internal solRad/d with astropy units. A file given in
  km/s is still written in m/s. EXOFASTv2 assumes m/s throughout.
- SED fluxes are lambda*F_lambda in erg/s/cm^2. Wavelengths are in microns.
- MIST tracks: Teff in K, radius in R_sun, age in Gyr.

## 3. The files, per component

Each subsection gives the EXOFASTv2 definition (file:line) and the EXOZIPPy quantity that
fills it. "At the data" means at the observed times, from the likelihood's own retained
node, never a second model (`components/instrument.md`, "The plotted model IS the
likelihood code").

### 3.1 Transit (`exofast_chi2v2.pro:1350-1357`, `plottran.pro:160-270`)

| file | columns (format) | EXOFASTv2 definition | EXOZIPPy source |
|---|---|---|---|
| `model.transit_%03d` | bjd, model (`%.8f %.6f`) | full model at the data: planets, F0, detrending, ramp, spline | the likelihood's mu at the data (detrend included, `_likelihood_mu_node`) **+ the GP conditional mean** (the GP stands where EXOFASTv2's spline does: inside the model, so the residuals are white) |
| `residuals.transit_%03d` | bjd, flux - model, err (`%.8f %.6f %.6f`) | data minus the line above, raw error | `self.flux - model`, `self.err` (raw, before `err_scale`/jitter, as EXOFASTv2 writes `transit.err`) |
| `detrendedmodel.transit_%03d.planet_%02d` | bjd, term (`%.10f %.10f`) | planet j's `exofast_tran - 1` at the data, exposure-smeared, before F0 and detrending | column j of `_eval_lc_data`'s `(N_obs, N_planets)` planet terms, i.e. planet j's step in the running model |
| `prettymodel.transit_%03d.planet_%02d` | time, term (`%.10f %.10f`) | the same, on a 1-per-minute grid | the same column from `Transit._lc_at_times(param_values, i, t)`, which smears over file i's own exposure as its data rows are |
| `spline.transit_%03d` | bjd, norm (`%.8f %.6f`) | the Kepler spline's multiplicative norm, only for `fitspline` files | `1 + (the spline columns' fitted term) / baseline`, only for `fitspline:` files (section 7) |

- The pretty grid is EXOFASTv2's, exactly: `lo = min(t_i) - 1`, `hi = max(t_i) + 1`,
  `n = ceil((hi - lo) * 1440)`, `t = lo + (hi - lo) * arange(n) / (n - 1)`.
- One planet file is written per (file, planet) for every planet the transit component
  models, as EXOFASTv2 does for every `fittran` planet, including planets with no event in
  that file (all zeros).
- EXOFASTv2's note 1 ("the detrended data can be derived by adding the residuals to the
  detrended model", i.e. `1 + sum_j detrendedmodel + residuals`) holds exactly for an
  additive baseline near 1 with no detrending and no GP. Otherwise it holds with the
  baseline in place of 1. Document that next to the file list rather than bending the
  numbers.

### 3.2 RV (`exofast_chi2v2.pro:1068-1092`, `plotrv.pro:95-215, 440-447`)

EXOZIPPy's RV model matrix has one column per orbit, Keplerian orbits first and then one
per Taylor orbit (`orbit/orbit.md`, "Taylor orbits"). A Taylor orbit is EXOFASTv2's
slope/quad (`/fitslope`, `/fitquad`), so "planets" below means the Keplerian columns and
"trend" means the Taylor columns (`system.orbit.is_taylor`).

| file | columns (format) | EXOFASTv2 definition | EXOZIPPy source |
|---|---|---|---|
| `model.telescope_%02d` | bjd, model (`%.8f %.6f`) | full model at the data: planets + gamma + slope/quad + detrending. (The header comment at `exofastv2.pro:1102` says gamma is subtracted. The code at `:1077-1088` includes it, and the code is what files contain.) | mu at the data (`_likelihood_mu_node`) + the GP conditional mean |
| `residuals.telescope_%02d` | bjd, rv - model, err, sigma (`%.8f %.6f %.6f %.6f`) | data minus the line above, raw error (upstream); `samsmods` writes `sqrt(err^2 + jitter variance)` instead | `rv - model`, raw `err`, then a **fourth column**, the likelihood's own total sigma at the point (decision 2) |
| `detrendedmodel.telescope_%02d` | bjd, model (`%.8f %.6f`) | planets only at the data, RM included; no gamma, trend or detrending | sum of the Keplerian columns of `_rv_matrix_data_node` at the data (an RM file's anomaly is on its own rows, as in the likelihood) |
| `trend.telescope_%02d` | bjd, trend (`%.8f %.6f`) | `gamma + slope (t - rvepoch) + quad (t - rvepoch)^2` at the data | gamma + the Taylor columns at the data (zero columns when there is no Taylor orbit: the file is then gamma, as in EXOFASTv2 with no slope) |
| `prettymodelrv.planet.%02d` | time, rv (`%.10f %.10f`) | planet i's curve on the pretty grid, no gamma/trend | planet i's orbit column (`planet.orbit_map[i]`) from `RVInstrument._rv_at_times` for the reference instrument (`_plot_layout`) |
| `prettymodelrv.trend` | time, planets, trend (`%.10f %.10f %.10f`) | col 2: `prettymodel`, the LAST planet's curve left over from the planet loop (`plotrv.pro:446`); col 3: slope/quad on the grid, no gamma | col 2: the sum of all Keplerian columns (decision 3; identical for a single planet); col 3: sum of the Taylor columns |

- The pretty grid is EXOFASTv2's, exactly (`plotrv.pro:73-99`): over all RV files,
  `roundto = 10**len(str(floor(max - min)))`, `bjd0 = floor(min / roundto) * roundto`,
  `lo = (min - bjd0) * 0.95 + bjd0`, `hi = (max - bjd0) * 1.05 + bjd0`, cadence = (shortest
  Keplerian period) / 100, `n = (hi - lo) / cadence`, `t = lo + (hi - lo) * arange(n) /
  (n - 1)`. Pin the integer conversion of `n` against IDL's `dindgen(nsteps)`, which
  truncates. This is not EXOZIPPy's own `_unphased_grid`, which serves the PDFs.
- A planet whose orbit has no RV signal (no RV file observes its star) gets no
  `prettymodelrv.planet` file, as EXOFASTv2 skips planets without `fitrv`.
- **The fourth residuals column** is `sigma_at_data` (section 4): the per-observation
  sigma the likelihood scored, `sqrt(err^2 + jitter_variance)` (`Instrument.total_sigma`),
  compiled from the likelihood's own node rather than recomputed. It is in m/s, written
  `%.6f`. EXOZIPPy's jitter variance is signed (`components/instrument.md`, "The signed
  jitter"), so where it is negative this column is SMALLER than the raw error. `samsmods`
  floors the jitter at zero and cannot do that. Do not clamp it: the column reports what
  the fit used. A reader of upstream's three columns (`readcol`, `np.loadtxt(usecols=...)`)
  is unaffected by the extra column.

### 3.3 Files EXOZIPPy adds (no EXOFASTv2 equivalent)

| file | columns (format) | content |
|---|---|---|
| `gp.transit_%03d`, `gp.telescope_%02d` | bjd, gp (`%.8f %.6f`) | the GP conditional mean at the data, in that file's units, for files with a GP. It is already inside `model.*`; this separates it, as `spline.transit_*` does in EXOFASTv2. |
| `<basename(prefix)>.index.txt` | text | dataset number -> name -> data file (section 2) |

### 3.4 SED (`sed/mistmultised.pro:460-490`, `exofast_chi2v2.pro:907-925`)

`<base>.sed.residuals.txt`, one row per photometric measurement (filter row), with
EXOFASTv2's comment line and format, verbatim:

```
# Filtername, Center wavelength (um), half bandpass (um), flux, error, modelflux, residuals (erg/s/cm^2), star indices
(a20,x,f0.6,x,f0.6,x,e0.6,x,e0.6,x,e0.6,x,e0.6,x,a)
```

- **Name and indices.** Filtername is the filter's name as the config spells it. Star
  indices are EXOFASTv2's blend string: the 0-based indices of the stars the row adds,
  comma-joined, then `-` and the stars it subtracts for a differential row (`0,1`, `1-0`,
  `-1`).
- **Wavelengths.** Center wavelength is the filter's effective wavelength. Half bandpass is
  half its effective width, from the same filter parameters the SED chart reads. If a
  filter has no effective width there, that RAISES, naming the filter.
- **Fluxes.** flux / error / modelflux / residuals are the observed and model
  lambda*F_lambda in erg/s/cm^2, from the SED component's own magnitude -> flux conversion
  (the one its chart draws in log10, and `sed_in_flux` on the summary figure). For a
  differential row both sides follow `_calc_obs_flux_from_obs_mag`, as the chart does.
  Error is `flux * ln(10) / 2.5 * mag_err`, with the raw error and no `errscale`, as
  EXOFASTv2 computes it.
- **Not identical numbers.** EXOFASTv2 converts with an AB zeropoint (`3631 Jy * 3e-9 /
  lambda_um`, Vega bands shifted by `vegaab`), while EXOZIPPy uses per-filter SVO
  zeropoints (`sed/make_bc.py`). The units and columns are identical; the fluxes differ at
  the zeropoint level. The writer does not imitate EXOFASTv2's zeropoints, because then the
  file would disagree with the fit it reports.

`<base>.atmosphere.%03d.txt`, one per star: wavelength (um) and that star's extincted model
lambda*F_lambda at Earth (erg/s/cm^2). EXOFASTv2 writes these in IDL free format, so use
`%.6f %.6e`. EXOFASTv2's grid is `0.1 + arange(24000) / 1000` um. EXOZIPPy writes the SED
component's own model spectrum on its native wavelengths (decision 5): the row count
differs from EXOFASTv2's, but the columns and units do not.

### 3.5 Evolutionary model (`massradius_mist.pro:313-341`)

`<base>.mist.%03d.track`, one per star with a MIST model, indexed per star even in a
single-star fit (decision 4; `samsmods` writes an unindexed `<base>.mist.track` beside the
PDFs). The columns are `teff rstar age
eep` (`%.5f` each), under the comment line `#teff, rstar, age, eep`, which is the 4-column
form upstream's `massradius_mist` writes. The rows are the whole track at the point's mass
and initial [Fe/H], `mist_grid.interpolate_track(grid, logmass, initfeh)`: every EEP whose
Teff is finite and whose radius is positive. That is EXOFASTv2's range (EEP 1-808,
`/allowold`, non-finite rows dropped). It does NOT apply the Kiel chart's EEP window or its
`here_be_dragons` cut, both of which are presentation choices for the PDF
(`MISTPlot._track_rows`). The writer reuses the interpolation, not `_track_rows`. A
non-MIST evolutionary model writes no track.

### 3.6 Everything else

Microlensing, astrometry, Doppler tomography and the pharmacokinetics components write
nothing: EXOFASTv2 has no text output for them. The hook (section 4) lets any of them add
files later.

## 4. Architecture: a component hook plus a generic writer

These formats are component knowledge: which quantity is "the model", what the planet
terms are, which units. So they live in the components ("All component-specific logic is
handled within the components", CLAUDE.md), and the writer is generic.

- **`Component.model_files(self, system, point) -> list[ModelFile]`**, an optional hook
  returning `[]` by default, like `summary_header`. `ModelFile` is a small dataclass in a
  new `src/exozippy/modelfile.py` (beside `chart.py`): `kind` (e.g. `"model.transit_000"`, or
  `"mist.000"` with `ext=".track"`), `columns` (a list of 1-D arrays, or strings for a
  text column), `formats` (one printf code per column), `comment` (a single header line,
  or None for `/nocomment`), `ext` (default `".txt"`) and `source` (the dataset's name
  and file path, for the index; None for a file that is not a dataset).
- **`outputs/model_files.py`: `write_model_files(system, point, prefix, tag)`** calls
  every active component's hook and writes `<base>.<kind><ext>`. It checks that the
  columns of one file have equal lengths, and that no two files share a name. Either
  failure RAISES, naming the component and the file. It then rewrites the index. It names
  no component; the `tests/test_summary_plot.py` yaml_key-literal scan is copied over to
  enforce that.
- **The likelihood's full mu and sigma at the data, in the base.** As built,
  `add_observation_likelihood` -- the one dispatcher every instrument's likelihood goes
  through -- retains the `mu` and `sigma` it is handed, and `Instrument.mu_at_data` /
  `sigma_at_data` compile them against `plot_params` LAZILY, on first use, so a run (or a
  GUI session) that never writes model files pays no extra compile; a rebuild resets the
  cache. (The plan said `compile_plotters` and `_likelihood_mu_node`; retaining at the
  dispatcher covers sigma too, and every child at once.) Those are the scored quantities,
  so the detrend comes from the one base mechanism (instrument.md, "One detrend
  mechanism") and no child rebuilds the forward map by hand. RV and transit both use them;
  `gp_mean_at_data` is added on top, as the plots do. `sigma_at_data` serves the RV
  residuals' fourth column. Both accessors are base-owned, for the reason the detrend
  mechanism is: a child cannot drift from its own likelihood.
- **Stale files.** Before writing a tag, remove that tag's previous files
  (`<base>.*` in `modelfiles/`), so a rerun with fewer datasets leaves no orphans.
  EXOFASTv2 does the same for the detrended models (`plottran.pro:162-163`). The glob is
  this prefix and tag only, so other fits that share the directory are untouched.

## 5. The `mcmc` point: the median draw

EXOFASTv2's `mcmc` files are its best-likelihood draw (`exofastv2.pro:1715-1719`).
EXOZIPPy writes them at `summary_plot.median_draw_point(system, posterior)`: one joint
draw, the valid draw nearest the posterior median in the sampled coordinates. That is the
point the summary figure draws, and it agrees with the medians in the results table. It is
not a vector of per-parameter medians, which need not be a point the posterior contains.
`outputs.md` ("The point is the median draw of the REPORTED posterior") records why the
highest-lp draw is wrong in EXOZIPPy. Each file therefore describes the reported fit, not
the luckiest draw.

## 6. Hooks in `run.py` and the console script

- **Start / polished** (`run.py`, around the `comp.plot(..., "_start")` and `"_polished"`
  loops, ~line 1167): `write_model_files(system, build_internal_start, prefix, "start")`,
  and the polished equivalent when `polished_start`.
- **Wrap-up** (`_wrap_up`, after `write_summary_plot`): a new `wrapup.stage("model
  files")` computes `median_draw_point(system, idata)` (idata is already trimmed,
  mode-labelled and distributed) and writes `mcmc`. Like every post-save stage it catches
  nothing (`run.md`, "A fit has three phases"). Add it to
  `tests/test_run_wrapup_guards.py`. `recompute_trace: false` reruns go through
  `_wrap_up`, so they rewrite the set.
- **Per mode** (`_emit_per_mode_outputs`): `mcmc_<mode_suffix(k)>` at
  `median_draw_point(system, _idata_for_mode(idata, k))`. Checked: `_idata_for_mode`
  returns a one-chain posterior of the mode's own draws, with no `mode` variable, which is
  the shape `median_draw_point` reads, so it shipped with the rest.
- **`exozippy-modelfiles <config.yaml>`** (`cli_modelfiles.py`, modeled on
  `cli_summary.py`), with `outputs.model_files.create_model_files(config_file)` modeled on
  `create_summary_plot`. It chdirs to the config's directory, builds the System
  (`prepare()` then `build_model()`), loads `<prefix>_trace.nc`, runs
  `check_trace_freshness`, `reported_posterior` and `median_draw_point`, then writes
  `mcmc`. `--start` adds the start set. There is no `--polished`: the polished point is
  the seed polish's result, i.e. the fit itself, and cannot be rebuilt from the config.
  Logging goes to the console only (the `<prefix>.log` reason in `cli_summary.py`).
  Registered in `[project.scripts]`. This gives existing fits their files without
  re-sampling, and its `mcmc` files are byte-identical to the ones the live wrap-up wrote
  (`tests/test_summary_plot.py`).

## 7. Known differences from EXOFASTv2 (documented, not bugs)

- `mcmc` is the median draw, not the best-likelihood draw (section 5).
- Dataset numbers follow the config's order (section 2; see the index file).
- SED fluxes use EXOZIPPy's zeropoints (section 3.4).
- The transit planet term is EXOZIPPy's additive step (`baseline + sum(terms) = model`).
  EXOFASTv2's is `exofast_tran - 1`, applied before multiplying by F0. They agree to
  `depth * (F0 - 1)`, about 1e-6 for a normalized light curve.
- `spline.transit_*` is `1 + term / baseline`, because EXOZIPPy's Kepler spline is additive
  detrend columns (instrument.md, "`fitspline`"). EXOFASTv2's is multiplicative. They
  agree to first order.
- The RV pretty curve is the reference instrument's (`_plot_layout`), which carries an RM
  anomaly only when every RV file has `rm:`. EXOFASTv2's pretty curve always includes RM,
  evaluated with `u1 = 0`.
- GP files are new; EXOFASTv2 has no GP.
- The polished start is tagged `polished`, not `amoeba`.
- `residuals.telescope_*` has a fourth column, the total sigma (section 3.2).
- `prettymodelrv.trend` column 2 is the sum of all planets, where EXOFASTv2 writes only the
  last planet's curve. The two are identical for a single planet.
- The `atmosphere.*` files are on the model's native wavelengths.

## 8. Decisions (settled 2026-10-14)

1. **The polished start is tagged `polished`, not `amoeba`.** The polish is not AMOEBA,
   and the PDFs already say `_polished`. The cost is that an EXOFASTv2 script globbing
   `*.amoeba.*` must be changed to `*.polished.*`. The `start` and `mcmc` tags match
   EXOFASTv2's.
2. **RV residuals: raw error, plus the jitter-inflated sigma as a fourth column.** Column 3
   is upstream's raw error, matching the transit files. Column 4 is what the fit actually
   used, which `samsmods` wrote in column 3's place (section 3.2: the likelihood's own
   sigma, signed jitter and all). Upstream's three-column readers keep working. The
   transit residuals stay at EXOFASTv2's three columns, which no branch extends.
3. **`prettymodelrv.trend` column 2 is the sum of all Keplerian orbits.** EXOFASTv2's
   column is the last planet's curve, a leftover from the planet loop. The sum is
   identical for a single planet and is evidently what the column means.
4. **MIST tracks are indexed per star**, `modelfiles/<base>.mist.%03d.track` even for a
   single star, after upstream's per-star `.mist.000.eps` naming. An unindexed name for
   star 0 only would break multi-star fits.
5. **Atmosphere files use the model's own wavelengths.** Interpolating onto EXOFASTv2's
   0.001-um grid would match the row count but invent resolution the model does not have.

## 9. Tests (as built)

`tests/test_model_files.py` (the writer), `tests/test_model_files_transit_rv.py`,
`tests/test_model_files_sed_mist.py`, and the end-to-end pair at the bottom of
`tests/test_summary_plot.py`, which reuse its finished kelt4 fit.

1. **Formats.** Every line of every transit, RV, SED and track file matches a regex built
   from EXOFASTv2's IDL format, with no header where EXOFASTv2 has none. The SED comment
   line and the track comment line are checked verbatim against the IDL sources. The
   `a20` field is right-justified as Fortran's and IDL's `A` edit descriptor writes on
   output. That is IDL's documented behaviour, NOT a comparison with a real
   `linkpars_current` file, because no such SED file was at hand. If one turns up, a
   one-line fixture would close that gap.
2. **Names.** The exact EXOFASTv2 names, padding included, for kelt4 (one TESS file, four
   RV files, one planet), and the index file. A duplicate file name, one dataset number
   naming two datasets, unequal column lengths and a format/column count mismatch each
   RAISE, and a raise leaves the previous set on disk intact.
3. **Numbers.** These are exact before formatting:
   - `model + residuals == data` for transit and RV.
   - `model.* == mu_at_data + gp`.
   - `baseline + detrendedmodel == model` for a light curve with no detrend.
   - `detrendedmodel.telescope + trend.telescope == model.telescope`.
   - The pretty transit curve at the data times equals `detrendedmodel`.
   - The pretty RV curve equals `_rv_at_times` for the reference instrument.
   - The spline file equals `1 + term / baseline`, at non-zero coefficients.
   - The SED's `log10(flux / model)` is the row's scored magnitude residual (with the
     differential side convention), its error is EXOFASTv2's.
   - The atmospheres are the chart's spectra.
   - The MIST track is `interpolate_track` at the Kiel node's mass and [Fe/H], past the
     chart's EEP window.

   `mu_at_data`/`sigma_at_data` reproduce PyMC's own logp of the observed RVs, and the RV
   fourth column is `sqrt(err^2 + jitter_variance)` and stays BELOW the raw error at a
   negative variance.
4. **Units.** The RV error column equals the RV file's own error column (m/s, with no
   constant in the test), and a km/s file is written in m/s. Times are the internal
   BJD_TDB the instruments hold; a `time_offset:` file is not separately tested.
5. **Grids.** The RV pretty grid is hand-checked against `plotrv.pro`'s arithmetic,
   including IDL's truncated point count over an untruncated `nsteps - 1`. The transit
   grid's count and endpoints are checked against `plottran.pro`'s.
6. **Component-agnostic writer.** The yaml_key-literal scan over `modelfile.py` and
   `outputs/model_files.py`, plus a made-up pharmacokinetics component's file.
7. **Run integration.** `tests/test_run_wrapup_guards.py` pins the wrap-up stage. The
   finished kelt4 fit has written its `start` and `mcmc` sets and the index, and
   `exozippy-modelfiles --start` rewrites the `mcmc` files byte-identically. The
   missing-trace message is not separately tested.
8. **Manual check, not CI.** Compare a fit's files with an EXOFASTv2 fit of the same
   system. Columns should agree to within the two fits' parameter differences, and the
   SED fluxes to within the zeropoint convention. A `samsmods` run has unpadded names and
   no subdirectory, so match by dataset, not by file name.

## 10. Docs

- `outputs/outputs.md`: a "Model files" section with the file table (sections 3.1-3.5 in
  brief), the reference branch, the point choice (pointing at the summary-figure rule
  rather than restating it), and section 7.
- `components/components.md`: the `model_files` hook in the component contract.
- `components/instrument.md`: `mu_at_data` and `sigma_at_data`, next to the retained
  nodes they compile.
- `run.md`: the new wrap-up stage and the start/polished sets.
- CLAUDE.md: the outputs line of the subsystem map names model files; README step 6
  documents `exozippy-modelfiles`.

## 11. Commit order

1. `ModelFile`, `outputs/model_files.py`, the `Component.model_files` default,
   `Instrument.mu_at_data` and `Instrument.sigma_at_data`, with tests 2 and 6. No component writes anything yet.
2. Transit and RV hooks, with tests 1 and 3-5 for them. This is the bulk of what users
   plot.
3. SED and MIST hooks, and the SED fixture line.
4. `run.py` hooks, `exozippy-modelfiles` and test 7.
5. Docs (section 10), or folded into each commit as it lands.
