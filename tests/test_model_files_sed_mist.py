"""
Tests for the SED and MIST model files (SED.model_files,
EvolutionaryModel.model_files): EXOFASTv2's sed.residuals, atmosphere and
mist track (docs/model_files_plan.md sections 3.4-3.5).

Two shipped builds: kelt4 transit+RV+SED (three stars, blended and
differential rows) and hat3 with MIST.

They follow AAA with Given/When/Then docstrings.
"""

import os
import re
from pathlib import Path

import numpy as np
import pytest
import yaml

from exozippy.components.evolutionarymodel import mist_grid
from exozippy.components.evolutionarymodel.mist_grid import OUTPUT_INDEX
from exozippy.components.evolutionarymodel.plot import MISTPlot
from exozippy.components.sed.sed import _exofast_star_indices
from exozippy.outputs.model_files import write_model_files
from exozippy.system import System

_EXAMPLES = Path(__file__).parent.parent / "examples"


def _build(directory, config_name):
    if not (_EXAMPLES / directory).is_dir():
        pytest.skip(f"{directory} example not present")
    cwd = os.getcwd()
    os.chdir(_EXAMPLES / directory)
    try:
        with open(config_name) as f:
            config = yaml.safe_load(f)
        with open(config.pop("parameter_file")) as f:
            params = yaml.safe_load(f)
        system = System(config, user_params=params)
        system.prepare()
        model = system.build_model()
        with model:
            point = system.get_internal_point(
                model, system.get_raw_start(model)
            )
    finally:
        os.chdir(cwd)
    return system, model, point


def _by_kind(files):
    return {f.kind: f for f in files}


def test_star_index_strings_are_exofasts():
    """
    Given a blend, a differential row and a row with only a negative side,
    When EXOFASTv2's star-index string is built,
    Then it is mistmultised.pro's: '0,1', '0-1,2', '-1'.
    """
    assert _exofast_star_indices([0, 1], []) == "0,1"
    assert _exofast_star_indices([0], [1, 2]) == "0-1,2"
    assert _exofast_star_indices([], [1]) == "-1"


@pytest.fixture(scope="module")
def kelt4():
    return _build("kelt4", "kelt4_rv+transit+sed.yaml")


@pytest.mark.slow
def test_sed_residuals_file_is_exofasts(kelt4, tmp_path):
    """
    Given the kelt4 build (three stars; Gaia rows of A alone, 2MASS/WISE
      rows of the A+B+C blend, a B-C and two A-(B+C) differential rows),
    When the model files are written,
    Then sed.residuals.txt starts with EXOFASTv2's comment line verbatim,
      has one line per measurement in the IDL format (a20, f0.6 x2, e0.6
      x4, a), and the star indices read 0 / 0,1,2 / 1-2 / 0-1,2.
    """
    system, _, point = kelt4
    sed = system.sed

    write_model_files(system, point, tmp_path / "kelt4", "mcmc")

    path = tmp_path / "modelfiles" / "kelt4.mcmc.sed.residuals.txt"
    lines = path.read_text().splitlines()
    assert lines[0] == (
        "# Filtername, Center wavelength (um), half bandpass (um), flux, "
        "error, modelflux, residuals (erg/s/cm^2), star indices"
    )
    rows = lines[1:]
    assert len(rows) == sed.nfilters
    e6 = r"-?\d\.\d{6}e[+-]\d{2}"
    regex = re.compile(rf"^.{{20}} \d+\.\d{{6}} \d+\.\d{{6}}( {e6}){{4}} \S+$")
    assert all(regex.match(r) for r in rows), rows[:3]
    names = [r[:20].strip() for r in rows]
    assert names == list(sed.filters)
    # (filter, stars) pairs: a filter may carry several rows (Keck Kp is
    # both the B-C and an A-(B+C) measurement).
    pairs = set(zip(names, (r.split()[-1] for r in rows)))
    assert {
        ("GAIA/GAIA3.G", "0"),
        ("2MASS/2MASS.J", "0,1,2"),
        ("Keck/NIRC2.Kp", "1-2"),
        ("Generic/Bessell.V", "0-1,2"),
    } <= pairs


@pytest.mark.slow
def test_sed_fluxes_carry_the_scored_magnitude_residual(kelt4):
    """
    Given the kelt4 build,
    When the sed.residuals columns are built,
    Then residual = flux - model, the error is flux * ln(10)/2.5 * the raw
      magnitude error, and log10(flux / model) is the row's magnitude
      residual against the model's combined prediction, -0.4 * side *
      (m_obs - m_pred) -- side +1 for a blend, -1 for a differential row,
      which EXOFASTv2 writes for its NEGATIVE (fainter) side.
    """
    system, _, point = kelt4
    sed = system.sed

    (f,) = [
        x for x in sed.model_files(system, point) if x.kind == "sed.residuals"
    ]

    _, _, _, flux, err, model, resid, _ = f.columns
    flux, err, model, resid = map(np.asarray, (flux, err, model, resid))
    np.testing.assert_allclose(resid, flux - model, rtol=1e-12)
    np.testing.assert_allclose(
        err, flux * np.log(10.0) / 2.5 * sed.err, rtol=1e-12
    )
    pred = sed._make_plot_obj(system, [point]).combined_pred_draws[0]
    side = np.where(np.any(sed.blend_matrix < 0, axis=1), -1.0, 1.0)
    np.testing.assert_allclose(
        np.log10(flux / model),
        -0.4 * side * (sed.mag - pred),
        rtol=0,
        atol=1e-10,
    )


@pytest.mark.slow
def test_atmosphere_files_are_the_charts_model_spectra(kelt4):
    """
    Given the kelt4 build (three stars),
    When the atmosphere files are built,
    Then there is one per star, atmosphere.000-002, on the model's own
      wavelengths, and each is the star's spectrum the SED chart draws.
    """
    system, _, point = kelt4
    sed = system.sed

    files = _by_kind(sed.model_files(system, point))

    (chart,) = sed.plot_data(system, point)
    spectra = {t.name: t for t in chart.traces if t.role == "model"}
    for s, name in enumerate(system.star.names):
        wave, lamflam = files[f"atmosphere.{s:03d}"].columns
        trace = spectra[f"Star {name}"]
        np.testing.assert_array_equal(wave, trace.x)
        ok = lamflam > 0
        np.testing.assert_allclose(
            np.log10(lamflam[ok]), np.asarray(trace.y)[ok], rtol=1e-12
        )
    assert f"atmosphere.{len(system.star.names):03d}" not in files


@pytest.fixture(scope="module")
def hat3_mist():
    return _build("hat3", "hat3_mist.yaml")


@pytest.mark.slow
def test_mist_track_is_the_whole_track_at_the_point(hat3_mist, tmp_path):
    """
    Given hat3 with a MIST model on star 0,
    When the model files are written,
    Then mist.000.track carries EXOFASTv2's comment line and teff, rstar,
      age, eep (%.5f) for every EEP the grid models along the track at the
      point's mass and initial [Fe/H] -- running past the Kiel chart's EEP
      window, which is the PDF's presentation, not the track.
    """
    system, _, point = hat3_mist
    em = system.evolutionarymodel

    write_model_files(system, point, tmp_path / "hat3", "start")

    path = tmp_path / "modelfiles" / "hat3.start.mist.000.track"
    lines = path.read_text().splitlines()
    assert lines[0] == "#teff, rstar, age, eep"
    f5 = r"-?\d+\.\d{5}"
    assert all(re.match(rf"^{f5} {f5} {f5} {f5}$", ln) for ln in lines[1:])
    got = np.loadtxt(path)
    kiel = np.atleast_2d(
        em._compiled_kiel(*em._point_to_plot_params(point, system))
    )
    track = mist_grid.interpolate_track(
        em._grids[0],
        float(kiel[0, MISTPlot.KIEL_INDEX["logmass"]]),
        float(kiel[0, MISTPlot.KIEL_INDEX["initfeh"]]),
    )
    eep = np.asarray(em._grids[0]["eep_pts"], dtype=float)
    keep = np.isfinite(track[:, OUTPUT_INDEX["teff_mist"]])
    np.testing.assert_allclose(got[:, 3], eep[keep])
    np.testing.assert_allclose(
        got[:, 0], track[keep, OUTPUT_INDEX["teff_mist"]], atol=5e-6
    )
    assert got[:, 3].max() > MISTPlot.KIEL_EEP_WINDOW[1]
