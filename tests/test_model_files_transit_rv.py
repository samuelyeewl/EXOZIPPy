"""
Tests for the transit and RV model files (Transit.model_files,
RVInstrument.model_files): EXOFASTv2's names, formats, grids and units, and
that every file's numbers are the likelihood's own (docs/model_files_plan.md
sections 3.1-3.3).

The numbers are checked on the ModelFile columns themselves (no formatting
round-off) and the formats on the files as written.  Two builds: the
shipped kelt4 transit+RV+SED example, and a small synthetic system that has
what kelt4 lacks -- a fitspline light curve, GPs on a light curve and on an
RV file, and an RV file given in km/s.

They follow AAA with Given/When/Then docstrings.
"""

import os
import re
from pathlib import Path
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
import yaml

from exozippy.components.rvinstrument.rvinstrument import RVInstrument
from exozippy.outputs.model_files import write_model_files
from exozippy.system import System

pytestmark = pytest.mark.slow

_KELT4_DIR = Path(__file__).parent.parent / "examples" / "kelt4"
_MS = float((u.solRad / u.d).to(u.m / u.s))


def _build(config, params=None):
    system = System(config, user_params=params)
    system.prepare()
    model = system.build_model()
    with model:
        point = system.get_internal_point(model, system.get_raw_start(model))
    return system, model, point


def _by_kind(files):
    return {f.kind: f for f in files}


@pytest.fixture(scope="module")
def kelt4():
    """kelt4 transit + RV + SED, built, at its start point."""
    if not _KELT4_DIR.is_dir():
        pytest.skip("kelt4 example not present")
    cwd = os.getcwd()
    os.chdir(_KELT4_DIR)
    try:
        with open("kelt4_rv+transit+sed.yaml") as f:
            config = yaml.safe_load(f)
        with open(config.pop("parameter_file")) as f:
            params = yaml.safe_load(f)
        return _build(config, params)
    finally:
        os.chdir(cwd)


# ---------------------------------------------------------------------------
# Names and formats
# ---------------------------------------------------------------------------

# One regex per file kind, from EXOFASTv2's IDL format: f0.8 -> 8 decimals.
_F = {n: rf"-?\d+\.\d{{{n}}}" for n in (6, 8, 10)}
_LINE = {
    "model": f"{_F[8]} {_F[6]}",
    "residuals.transit": f"{_F[8]} {_F[6]} {_F[6]}",
    "residuals.telescope": f"{_F[8]} {_F[6]} {_F[6]} {_F[6]}",
    "detrendedmodel.transit": f"{_F[10]} {_F[10]}",
    "prettymodel": f"{_F[10]} {_F[10]}",
    "detrendedmodel.telescope": f"{_F[8]} {_F[6]}",
    "trend.telescope": f"{_F[8]} {_F[6]}",
    "prettymodelrv.planet": f"{_F[10]} {_F[10]}",
    "prettymodelrv.trend": f"{_F[10]} {_F[10]} {_F[10]}",
}


def _line_regex(kind):
    for key in sorted(_LINE, key=len, reverse=True):
        if kind.startswith(key):
            return re.compile(f"^{_LINE[key]}$")
    raise AssertionError(f"no format expectation for {kind!r}")


def test_kelt4_files_have_exofast_names_and_formats(kelt4, tmp_path):
    """
    Given the kelt4 build (one TESS file, four RV files, one planet),
    When its transit and RV model files are written under the start tag,
    Then the names are EXOFASTv2's, padding included (transit_%03d,
      telescope_%02d, planet_%02d), and every line of every file matches
      its IDL format with no header -- EXOFASTv2's /nocomment.
    """
    system, _, point = kelt4

    written = write_model_files(system, point, tmp_path / "kelt4", "start")

    names = {p.name for p in written}
    expected = {"kelt4.index.txt"}
    expected |= {
        f"kelt4.start.{k}.transit_000.txt" for k in ("model", "residuals")
    }
    expected |= {
        f"kelt4.start.{k}.transit_000.planet_00.txt"
        for k in ("detrendedmodel", "prettymodel")
    }
    expected |= {
        f"kelt4.start.{k}.telescope_{i:02d}.txt"
        for k in ("model", "residuals", "detrendedmodel", "trend")
        for i in range(4)
    }
    expected |= {
        "kelt4.start.prettymodelrv.planet.00.txt",
        "kelt4.start.prettymodelrv.trend.txt",
    }
    assert expected <= names
    for p in written:
        if p.name == "kelt4.index.txt" or not (
            ".transit_" in p.name or "telescope" in p.name or "rv." in p.name
        ):
            continue
        kind = p.name[len("kelt4.start.") : -len(".txt")]
        regex = _line_regex(kind)
        lines = p.read_text().splitlines()
        assert lines, p.name
        bad = [ln for ln in lines if not regex.match(ln)]
        assert not bad, (p.name, bad[:3])
    index = (
        (tmp_path / "modelfiles" / "kelt4.index.txt").read_text().splitlines()
    )
    assert index[0].startswith("transit_000  TESS")
    assert [ln.split()[0] for ln in index[1:]] == [
        f"telescope_{i:02d}" for i in range(4)
    ]


# ---------------------------------------------------------------------------
# RV numbers
# ---------------------------------------------------------------------------


def test_rv_files_add_up_to_the_scored_model_in_metres_per_second(kelt4):
    """
    Given the kelt4 build (no RV detrend columns, no RV GP, no Taylor
      orbit),
    When its RV model files are built,
    Then model = detrendedmodel + trend, model + residual = the data, the
      model is the likelihood's mu, and the error column is each file's own
      raw error column as read from disk (m/s, the files' unit) -- a unit
      check that carries no conversion constant of its own.
    """
    system, _, point = kelt4
    rv = system.rvinstrument
    assert rv.total_detrend_cols == 0 and not rv.has_gp

    files = _by_kind(rv.model_files(system, point))

    mu = rv.mu_at_data(system, point)
    for i in range(rv.n_elements):
        tag = f"telescope_{i:02d}"
        rows = rv.rows(i)
        t, model = files[f"model.{tag}"].columns
        _, resid, err, _ = files[f"residuals.{tag}"].columns
        _, planets = files[f"detrendedmodel.{tag}"].columns
        _, trend = files[f"trend.{tag}"].columns
        np.testing.assert_allclose(model, planets + trend, rtol=1e-12)
        np.testing.assert_allclose(model, mu[rows] * _MS, rtol=1e-12)
        np.testing.assert_allclose(model + resid, rv.rv[rows] * _MS)
        raw = np.loadtxt(_KELT4_DIR / rv.files[i])
        raw = raw[np.argsort(raw[:, 0])]
        np.testing.assert_allclose(t, raw[:, 0])
        np.testing.assert_allclose(err, raw[:, 2], rtol=1e-12)
        np.testing.assert_allclose(model + resid, raw[:, 1], rtol=1e-10)


def test_rv_fourth_column_is_the_scored_sigma_and_is_never_clamped(kelt4):
    """
    Given the kelt4 build at a point whose first RV jitter variance is
      NEGATIVE (the signed jitter: the quoted errors are too large),
    When its RV residuals are built,
    Then column 4 is sqrt(err^2 + jitter variance) -- smaller than the raw
      error in column 3, as the likelihood scored it, not floored at it.
    """
    system, _, point = kelt4
    rv = system.rvinstrument
    jv = np.array(rv._point_value(point, rv.jitter_variance), dtype=float)
    err0 = rv.err[rv.rows(0)]
    jv[0] = -0.5 * float(np.min(err0)) ** 2
    shifted = {**point, rv.jitter_variance.label: jv}

    files = _by_kind(rv.model_files(system, shifted))

    _, _, err, sigma = files["residuals.telescope_00"].columns
    np.testing.assert_allclose(
        sigma, np.sqrt(err0**2 + jv[0]) * _MS, rtol=1e-12
    )
    assert np.all(sigma < err)


def test_rv_pretty_grid_is_exofasts_hand_checked():
    """
    Given RV epochs spanning 300 days from BJD 2456000 and a shortest
      period of 200 days,
    When EXOFASTv2's pretty grid is built,
    Then by plotrv.pro's arithmetic: 3 digits -> round to 1000, bjd0 =
      2456000, the span widens to [2456000, 2456315], the cadence is 2 d so
      nsteps = 157.5 -> 157 points, spaced 315 / 156.5 (IDL divides by the
      untruncated nsteps - 1), so the last point falls short of 2456315.
    """
    stub = SimpleNamespace(time=np.array([2456000.0, 2456150.0, 2456300.0]))

    t = RVInstrument._exofast_pretty_rv_grid(stub, 200.0)

    assert t.size == 157
    assert t[0] == 2456000.0
    np.testing.assert_allclose(np.diff(t), 315.0 / 156.5)
    np.testing.assert_allclose(t[-1], 2456000.0 + 315.0 * 156 / 156.5)


def test_rv_pretty_planet_curve_is_the_reference_instruments_model(kelt4):
    """
    Given the kelt4 build (one planet, no Taylor orbit),
    When its pretty RV files are built,
    Then the planet curve is the reference instrument's plotted model on
      that grid, and the trend file's planets column is that same curve and
      its trend column zero.
    """
    system, _, point = kelt4
    rv = system.rvinstrument

    files = _by_kind(rv.model_files(system, point))

    t, curve = files["prettymodelrv.planet.00"].columns
    params = rv._point_to_plot_params(point, system)
    full, _ = rv._rv_at_times(params, rv._rv_layout[0], t)
    np.testing.assert_allclose(curve, full * _MS, rtol=1e-12)
    _, planets, trend = files["prettymodelrv.trend"].columns
    np.testing.assert_array_equal(planets, curve)
    assert np.all(trend == 0.0)


# ---------------------------------------------------------------------------
# Transit numbers
# ---------------------------------------------------------------------------


def test_transit_files_add_up_to_the_scored_model(kelt4):
    """
    Given the kelt4 build (one TESS file, no detrend columns, no GP),
    When its transit model files are built,
    Then model = the likelihood's mu = baseline + the planet's term,
      model + residual = the flux, the error column is the raw error, and
      the pretty curve evaluated at the data times is the detrendedmodel
      -- the two files describe one curve.
    """
    system, _, point = kelt4
    tr = system.transit
    assert tr.total_detrend_cols == 0 and not tr.has_gp

    files = _by_kind(tr.model_files(system, point))

    rows = tr.rows(0)
    t, model = files["model.transit_000"].columns
    _, resid, err = files["residuals.transit_000"].columns
    _, term = files["detrendedmodel.transit_000.planet_00"].columns
    baseline = tr._point_value(point, tr.baseline, 0)
    np.testing.assert_allclose(
        model, tr.mu_at_data(system, point)[rows], rtol=1e-12
    )
    np.testing.assert_allclose(model, baseline + term, rtol=1e-12)
    np.testing.assert_allclose(model + resid, tr.flux[rows], rtol=1e-12)
    np.testing.assert_array_equal(err, tr.err[rows])
    params = tr._point_to_plot_params(point, system)
    _, at_data = tr._lc_at_times(params, 0, t)
    np.testing.assert_allclose(at_data[:, 0], term, atol=1e-12)


def test_transit_pretty_grid_is_one_per_minute_a_day_past_the_data(kelt4):
    """
    Given the kelt4 TESS file,
    When its pretty model is built,
    Then its grid runs from a day before the first point to a day after the
      last, with ceil(span * 1440) points -- plottran.pro's grid.
    """
    system, _, point = kelt4
    tr = system.transit

    files = _by_kind(tr.model_files(system, point))

    t_data = tr.time[tr.rows(0)]
    t, _ = files["prettymodel.transit_000.planet_00"].columns
    lo, hi = t_data.min() - 1.0, t_data.max() + 1.0
    assert t.size == int(np.ceil((hi - lo) * 1440.0))
    assert t[0] == lo
    np.testing.assert_allclose(t[-1], hi, rtol=0, atol=1e-9)


# ---------------------------------------------------------------------------
# Spline, GP and a km/s RV file (synthetic)
# ---------------------------------------------------------------------------

_TC = 2459000.5
_PERIOD = 3.0


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    """Two light curves -- A with fitspline and an SHO GP, B plain -- and
    one RV file in km/s with an SHO GP, around one planet."""
    d = tmp_path_factory.mktemp("synthetic")
    rng = np.random.default_rng(1)
    t = np.linspace(_TC - 0.4, _TC + 0.4, 300)
    for name, wobble in (("a", 2e-3), ("b", 0.0)):
        flux = 1.0 + wobble * np.sin(6.0 * (t - _TC))
        flux += rng.normal(0, 5e-4, t.size)
        np.savetxt(
            d / f"{name}.dat",
            np.column_stack([t, flux, np.full(t.size, 5e-4)]),
        )
    t_rv = _TC + np.sort(rng.uniform(-30, 30, 25))
    rv_kms = 0.05 * np.sin(2 * np.pi * (t_rv - _TC) / _PERIOD)
    rv_kms += rng.normal(0, 0.005, t_rv.size)
    np.savetxt(
        d / "c.rv", np.column_stack([t_rv, rv_kms, np.full(t_rv.size, 0.005)])
    )
    config = {
        "star": [{"name": "A", "mist": False}],
        "planet": [{"name": "b"}],
        "orbit": [{"name": "b"}],
        "band": [{"name": "TESS", "filter": "TESS", "ld_law": "quadratic"}],
        "transit": [
            {
                "name": "A",
                "file": str(d / "a.dat"),
                "band": "TESS",
                "fitspline": True,
                "gp": "sho",
            },
            {"name": "B", "file": str(d / "b.dat"), "band": "TESS"},
        ],
        "rvinstrument": [
            {"name": "C", "file": str(d / "c.rv"), "unit": "km/s", "gp": "sho"}
        ],
    }
    params = {
        "star.0.radius": {"initval": 1.0, "sigma": 0.0},
        "star.0.logmass": {"initval": 0.0, "sigma": 0.0},
        "orbit.0.logP": {"initval": float(np.log10(_PERIOD)), "sigma": 0.0},
        "orbit.0.tc": {"initval": _TC, "sigma": 0.0},
    }
    system, model, point = _build(config, params)
    # A non-zero spline: at the start every detrend coefficient is 0.
    tr = system.transit
    coeffs = np.array(tr._point_value(point, tr.detrend_coeffs))
    coeffs = rng.normal(0.0, 1e-3, coeffs.size)
    point = {**point, tr.detrend_coeffs.label: coeffs}
    return system, model, point, d


def test_spline_file_is_one_plus_the_spline_term_over_the_baseline(
    synthetic,
):
    """
    Given light curve A with fitspline (its only detrend columns are the
      spline's) at non-zero coefficients, and plain light curve B,
    When the transit model files are built,
    Then A gets spline.transit_000 = 1 + (its detrend term) / baseline,
      EXOFASTv2's multiplicative norm to first order, and B gets none.
    """
    system, _, point, _ = synthetic
    tr = system.transit
    assert tr.n_spline_cols[0] == tr.n_detrend_per_inst[0] > 0

    files = _by_kind(tr.model_files(system, point))

    rows = tr.rows(0)
    _, spline_norm = files["spline.transit_000"].columns
    baseline = tr._point_value(point, tr.baseline, 0)
    term = tr.detrend_at_data(point)[rows]
    assert np.ptp(term) > 0
    np.testing.assert_allclose(spline_norm, 1.0 + term / baseline, rtol=1e-12)
    assert "spline.transit_001" not in files


def test_gp_is_inside_the_model_and_written_on_its_own(synthetic):
    """
    Given light curve A and RV file C with SHO GPs, and B without,
    When the model files are built,
    Then A's and C's model is the likelihood's mu plus the GP conditional
      mean, their gp file is that mean, and B gets no gp file.
    """
    system, _, point, _ = synthetic
    tr, rv = system.transit, system.rvinstrument

    lc = _by_kind(tr.model_files(system, point))
    rvf = _by_kind(rv.model_files(system, point))

    gp = tr.gp_mean_at_data(system, point)
    rows = tr.rows(0)
    assert np.ptp(gp[rows]) > 0
    _, model = lc["model.transit_000"].columns
    _, gp_col = lc["gp.transit_000"].columns
    np.testing.assert_array_equal(gp_col, gp[rows])
    np.testing.assert_allclose(
        model, tr.mu_at_data(system, point)[rows] + gp[rows], rtol=1e-12
    )
    assert "gp.transit_001" not in lc
    gp_rv = rv.gp_mean_at_data(system, point)
    _, rv_model = rvf["model.telescope_00"].columns
    _, rv_gp = rvf["gp.telescope_00"].columns
    np.testing.assert_allclose(rv_gp, gp_rv * _MS, rtol=1e-12)
    np.testing.assert_allclose(
        rv_model, (rv.mu_at_data(system, point) + gp_rv) * _MS, rtol=1e-12
    )


def test_an_rv_file_in_kms_is_written_in_metres_per_second(synthetic):
    """
    Given RV file C, given in km/s,
    When its residuals are built,
    Then the error column is the file's error column converted to m/s.
    """
    system, _, point, d = synthetic
    rv = system.rvinstrument

    files = _by_kind(rv.model_files(system, point))

    raw = np.loadtxt(d / "c.rv")
    _, _, err, _ = files["residuals.telescope_00"].columns
    kms = float((u.km / u.s).to(u.m / u.s))
    np.testing.assert_allclose(err, raw[:, 2] * kms, rtol=1e-12)
