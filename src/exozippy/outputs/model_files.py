"""
Write the EXOFASTv2-compatible model files of a fit at one point.

Every active component's ``model_files(system, point)`` hook returns the
files it owns (``exozippy.modelfile.ModelFile``); this module only places
them, exactly where EXOFASTv2 does:

    <dirname(prefix)>/modelfiles/<basename(prefix)>.<tag>.<kind><ext>

e.g. ``fitresults/modelfiles/toi2876.mcmc.model.transit_000.txt``.  It also
writes ``<basename(prefix)>.index.txt`` there, which maps each numbered
dataset (``transit_000``) to its name and data file -- EXOFASTv2's files have
no header, and EXOZIPPy numbers datasets in config order, not EXOFASTv2's
glob order.

Like the summary figure, this module names no component
(``tests/test_model_files.py`` scans it for component yaml_keys): a component
from any field joins by returning files from its hook.  The file definitions,
the EXOFASTv2 branch they follow and the decisions behind them are in
``docs/model_files_plan.md``.
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

#: The subdirectory EXOFASTv2 writes its model files to (exofastv2.pro,
#: ``modeldirname``).
MODEL_DIR = "modelfiles"


def model_dir(prefix):
    """``<dirname(prefix)>/modelfiles`` -- where a fit's model files go."""
    return Path(prefix).parent / MODEL_DIR


def write_model_files(system, point, prefix, tag):
    """Write every component's model files at ``point`` under ``tag``.

    ``point`` is an internal-unit point dict, as ``plot_data`` takes;
    ``tag`` is the EXOFASTv2-style stage name (``start``, ``polished``,
    ``mcmc``, ``mcmc_<mode>``).  The tag's previous files are removed first,
    so a rerun with fewer datasets leaves no orphans (EXOFASTv2 does the
    same for its detrended models, plottran.pro).  Returns the written
    paths, the index file last.

    RAISES when two components' files share a name, or when a component's
    file is malformed (``ModelFile.lines``): either is a bug upstream.
    """
    stem = Path(prefix).name
    base = f"{stem}.{tag}"

    # Every file is rendered before anything on disk changes, so a
    # malformed file raises with the previous set still intact.
    texts = {}
    owners = {}
    datasets = {}
    for key, comp in system.active_components.items():
        for mf in comp.model_files(system, point):
            name = f"{base}.{mf.kind}{mf.ext}"
            if name in owners:
                raise ValueError(
                    f"model file {name!r} is written by both "
                    f"{owners[name]!r} and {key!r}."
                )
            owners[name] = key
            texts[name] = "\n".join(mf.lines()) + "\n"
            if mf.dataset is not None:
                ident, description = mf.dataset
                seen = datasets.setdefault(ident, description)
                if seen != description:
                    raise ValueError(
                        f"model files number two datasets {ident!r}: "
                        f"{seen!r} and {description!r}."
                    )

    if not texts:
        logger.info(
            "No component of this fit writes model files; nothing written "
            "for %s.",
            tag,
        )
        return []

    outdir = model_dir(prefix)
    outdir.mkdir(parents=True, exist_ok=True)
    for old in outdir.glob(f"{base}.*"):
        old.unlink()

    written = []
    for name, text in texts.items():
        path = outdir / name
        path.write_text(text)
        written.append(path)

    index = outdir / f"{stem}.index.txt"
    index.write_text(
        "".join(f"{ident}  {desc}\n" for ident, desc in datasets.items())
    )
    written.append(index)
    logger.info(
        "Wrote %d model files (%s) to %s", len(texts), tag, outdir.resolve()
    )
    return written


def create_model_files(config_file, *, start=False):
    """Write a finished fit's model files from its config and saved trace.

    ``config_file`` is the system YAML the fit ran with: it rebuilds the
    System and its ``prefix:`` locates ``<prefix>_trace.nc``, which must
    belong to the model the config builds (``check_trace_freshness``).
    Relative paths resolve against the config's own directory, wherever
    this is called from -- as ``summary_plot.create_summary_plot`` does.

    Writes the ``mcmc`` set at the median draw of the reported posterior
    (``summary_plot.reported_posterior`` + ``median_draw_point``, exactly
    what the live fit's wrap-up writes), and with ``start=True`` also the
    ``start`` set at the build start.  The ``polished`` set needs the seed
    polish, i.e. the fit itself, and is not reproducible from here.
    Returns the written paths (absolute).
    """
    import contextlib

    import arviz as az

    from ..system import System
    from ..trace_meta import check_trace_freshness
    from ..yamlio import load_system_config
    from .summary_plot import (
        DEFAULT_PREFIX,
        median_draw_point,
        reported_posterior,
    )

    config_path = Path(config_file).expanduser().resolve()
    with contextlib.chdir(config_path.parent):
        config = load_system_config(config_path.name)
        prefix = Path(config.get("prefix", DEFAULT_PREFIX))
        trace_path = Path(f"{prefix}_trace.nc")
        if not trace_path.exists():
            raise FileNotFoundError(
                f"No saved trace at {trace_path.resolve()}.  The mcmc model "
                f"files are written from a finished fit; run `exozippy "
                f"{config_path.name}` first."
            )

        system = System(config)
        system.prepare()
        model = system.build_model()

        written = []
        if start:
            build_start = system.get_internal_point(
                model, system.get_raw_start(model)
            )
            written += write_model_files(system, build_start, prefix, "start")

        idata = az.from_netcdf(str(trace_path))
        check_trace_freshness(idata, system, trace_path)
        posterior = reported_posterior(system, model, idata)
        point, _, _ = median_draw_point(system, posterior)
        written += write_model_files(system, point, prefix, "mcmc")
        return [p.resolve() for p in written]
