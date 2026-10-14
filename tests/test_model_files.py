"""
Tests for the EXOFASTv2-compatible model files (exozippy/modelfile.py,
outputs/model_files.py, Component.model_files, Instrument.mu_at_data /
sigma_at_data).

The writer only places files; what a file HOLDS is its component's, so the
per-component numbers are tested where the files are built.  Here: the
rendering of one file (IDL formats translated literally), where the writer
puts files and what it refuses, that it names no component, and that the
likelihood's own mean and scale come back out of a built model exactly.
See docs/model_files_plan.md.

They follow AAA with Given/When/Then docstrings.
"""

import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from exozippy.modelfile import E6, F6, F8, ModelFile, a_field
from exozippy.outputs import model_files as mf_out

_KELT4_DIR = Path(__file__).parent.parent / "examples" / "kelt4"


# ---------------------------------------------------------------------------
# One file
# ---------------------------------------------------------------------------


def test_a_file_renders_exofast_formats_with_no_header_by_default():
    """
    Given a residuals-style file (bjd %.8f, value %.6f, error %.6f) and no
      comment,
    When its lines are rendered,
    Then each row is the three numbers in exactly those formats, one space
      apart, with no header line -- EXOFASTv2's /nocomment.
    """
    f = ModelFile(
        kind="residuals.transit_000",
        columns=[
            [2458491.765098, 2458491.785932],
            [2.51e-4, -1.3e-4],
            [7.84e-4, 7.82e-4],
        ],
        formats=[F8, F6, F6],
    )

    lines = f.lines()

    assert lines == [
        "2458491.76509800 0.000251 0.000784",
        "2458491.78593200 -0.000130 0.000782",
    ]


def test_a_text_column_takes_idl_a_formats_and_a_comment_is_verbatim():
    """
    Given an SED-style file: a 20-wide filter name (IDL a20), a wavelength,
      an exponential flux and a free-width text column, under a comment,
    When it is rendered,
    Then the comment is the first line unchanged, the name is right-
      justified in 20 characters (Fortran/IDL A output), the flux is %.6e.
    """
    f = ModelFile(
        kind="sed.residuals",
        columns=[["Bessell_B"], [0.436], [1.2219852e-10], ["0"]],
        formats=[a_field(20), F6, E6, "%s"],
        comment="# Filtername, flux",
    )

    lines = f.lines()

    assert lines == [
        "# Filtername, flux",
        "           Bessell_B 0.436000 1.221985e-10 0",
    ]


def test_a_malformed_file_raises_rather_than_padding():
    """
    Given a file whose columns differ in length, and one whose format count
      differs from its column count,
    When either is rendered,
    Then it raises naming the file: both are a component bug.
    """
    ragged = ModelFile("model.x_00", [[1.0, 2.0], [1.0]], [F8, F6])
    short = ModelFile("model.x_01", [[1.0], [1.0]], [F8])

    with pytest.raises(ValueError, match="model.x_00"):
        ragged.lines()
    with pytest.raises(ValueError, match="model.x_01"):
        short.lines()


# ---------------------------------------------------------------------------
# The writer
# ---------------------------------------------------------------------------


def _system(**hooks):
    """A stub system whose components are just model_files hooks."""
    return SimpleNamespace(
        active_components={
            key: SimpleNamespace(model_files=fn) for key, fn in hooks.items()
        }
    )


def _one(kind, dataset=None):
    return ModelFile(kind, [[1.0], [2.0]], [F8, F6], dataset=dataset)


def test_files_go_to_modelfiles_beside_the_prefix_with_an_index(tmp_path):
    """
    Given a prefix fitresults/toi2876 and a component writing one numbered
      dataset's model,
    When the mcmc set is written,
    Then it lands at fitresults/modelfiles/toi2876.mcmc.model.transit_000.txt,
      and fitresults/modelfiles/toi2876.index.txt maps transit_000 to its
      description.
    """
    prefix = tmp_path / "fitresults" / "toi2876"
    system = _system(
        transit=lambda s, p: [
            _one("model.transit_000", ("transit_000", "TESS  a.dat"))
        ]
    )

    written = mf_out.write_model_files(system, {}, prefix, "mcmc")

    outdir = tmp_path / "fitresults" / "modelfiles"
    assert written == [
        outdir / "toi2876.mcmc.model.transit_000.txt",
        outdir / "toi2876.index.txt",
    ]
    assert written[0].read_text() == "1.00000000 2.000000\n"
    assert written[1].read_text() == "transit_000  TESS  a.dat\n"


def test_a_rewrite_removes_only_that_tags_previous_files(tmp_path):
    """
    Given a start set and an mcmc set already on disk, the mcmc set holding
      a file the next run will not write, and another fit's file in the
      same directory,
    When the mcmc set is written again,
    Then the orphaned mcmc file is gone, and the start set and the other
      fit's file are untouched.
    """
    prefix = tmp_path / "fit"
    outdir = tmp_path / "modelfiles"
    outdir.mkdir()
    keep = [
        outdir / "fit.start.model.transit_000.txt",
        outdir / "other.mcmc.model.transit_000.txt",
    ]
    orphan = outdir / "fit.mcmc.model.transit_001.txt"
    for p in keep + [orphan]:
        p.write_text("x\n")
    system = _system(transit=lambda s, p: [_one("model.transit_000")])

    mf_out.write_model_files(system, {}, prefix, "mcmc")

    assert not orphan.exists()
    assert all(p.read_text() == "x\n" for p in keep)
    assert (outdir / "fit.mcmc.model.transit_000.txt").exists()


def test_two_components_writing_one_name_raise_and_nothing_changes(tmp_path):
    """
    Given two components that both return model.transit_000, and a
      previous mcmc set on disk,
    When the mcmc set is written,
    Then it raises naming both components, and the previous set is intact:
      every file is rendered before anything on disk changes.
    """
    prefix = tmp_path / "fit"
    outdir = tmp_path / "modelfiles"
    outdir.mkdir()
    old = outdir / "fit.mcmc.model.transit_000.txt"
    old.write_text("old\n")
    system = _system(
        a=lambda s, p: [_one("model.transit_000")],
        b=lambda s, p: [_one("model.transit_000")],
    )

    with pytest.raises(ValueError, match="'a' and 'b'"):
        mf_out.write_model_files(system, {}, prefix, "mcmc")

    assert old.read_text() == "old\n"


def test_one_dataset_number_naming_two_datasets_raises(tmp_path):
    """
    Given two files that number different datasets transit_000,
    When they are written,
    Then it raises: the index would otherwise point one number at two files.
    """
    system = _system(
        transit=lambda s, p: [
            _one("model.transit_000", ("transit_000", "A  a.dat")),
            _one("residuals.transit_000", ("transit_000", "B  b.dat")),
        ]
    )

    with pytest.raises(ValueError, match="transit_000"):
        mf_out.write_model_files(system, {}, tmp_path / "fit", "mcmc")


def test_a_fit_with_no_model_files_writes_nothing(tmp_path):
    """
    Given a fit none of whose components writes model files (a
      microlensing fit),
    When the set is written,
    Then nothing is written and no modelfiles/ directory is created.
    """
    system = _system(mulens=lambda s, p: [])

    written = mf_out.write_model_files(system, {}, tmp_path / "fit", "mcmc")

    assert written == []
    assert not (tmp_path / "modelfiles").exists()


def test_the_writer_names_no_component():
    """
    Given every component the factory discovers,
    When the writer's source is read,
    Then no string literal in it is a component's yaml_key: what a file
      holds is its component's, and the writer only places it.
    """
    import ast
    import inspect

    from exozippy import modelfile
    from exozippy.components.factory import discover_components

    keys = set(discover_components())
    literals = set()
    for module in (mf_out, modelfile):
        tree = ast.parse(inspect.getsource(module))
        literals |= {
            node.value
            for node in ast.walk(tree)
            if isinstance(node, ast.Constant) and isinstance(node.value, str)
        }

    assert keys and not (literals & keys), sorted(literals & keys)


def test_a_component_from_another_field_joins_by_returning_files(tmp_path):
    """
    Given a system whose only component is not an astronomy one (yaml_key
      "assay"), returning one concentration file,
    When the model files are written,
    Then that file is written like any other.
    """
    system = _system(
        assay=lambda s, p: [
            ModelFile(
                "model.subject_00",
                [[0.5, 1.0], [3.2, 2.9]],
                ["%.2f", "%.3f"],
                dataset=("subject_00", "patient 1  theoph.csv"),
            )
        ]
    )

    written = mf_out.write_model_files(system, {}, tmp_path / "pk", "mcmc")

    assert written[0].name == "pk.mcmc.model.subject_00.txt"
    assert written[0].read_text() == "0.50 3.200\n1.00 2.900\n"


# ---------------------------------------------------------------------------
# The likelihood's own mean and scale
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def rvonly_built():
    """kelt4 RV-only system, prepared and built, with a start point."""
    if not _KELT4_DIR.is_dir():
        pytest.skip("kelt4 example not present")
    cwd = os.getcwd()
    os.chdir(_KELT4_DIR)
    try:
        with open("kelt4_rvonly.yaml") as f:
            config = yaml.safe_load(f)
        from exozippy.system import System

        system = System(config)
        system.prepare()
        model = system.build_model()
        with model:
            point = system.get_internal_point(
                model, system.get_raw_start(model)
            )
    finally:
        os.chdir(cwd)
    return system, model, point


@pytest.mark.slow
def test_mu_and_sigma_at_data_score_exactly_what_the_likelihood_scores(
    rvonly_built,
):
    """
    Given the kelt4 RV-only build at its start point,
    When the Gaussian log-likelihood of the RVs is computed from
      mu_at_data and sigma_at_data,
    Then it equals the observed variable's own logp in the PyMC model --
      the two accessors return the scored mean and scale, not a rebuild of
      them -- and sigma is the quoted error in quadrature with the jitter.
    """
    from scipy.stats import norm

    system, model, point = rvonly_built
    rv = system.rvinstrument

    mu = rv.mu_at_data(system, point)
    sigma = rv.sigma_at_data(system, point)

    obs = model[f"{rv.prefix}.model"]
    logp_fn = model.compile_logp(vars=[obs], sum=True)
    with model:
        raw = system.get_raw_start(model)
    expected = float(logp_fn(raw))
    got = float(np.sum(norm.logpdf(rv.rv, loc=mu, scale=sigma)))
    # ~1e-10 apart, not bit for bit: the point comes out of one compiled
    # graph (raw -> deterministics, get_internal_point) and mu out of a
    # second, while PyMC's logp is one fused graph, so they round
    # differently.  A wrong term -- a dropped jitter or detrend -- moves this
    # logp by O(1).
    np.testing.assert_allclose(got, expected, rtol=1e-8)

    jv = rv._point_value(point, rv.jitter_variance)[rv.inst_map]
    np.testing.assert_allclose(sigma, np.sqrt(rv.err**2 + jv), rtol=1e-12)
