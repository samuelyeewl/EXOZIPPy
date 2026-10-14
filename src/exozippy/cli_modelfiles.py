"""Console entry point: exozippy-modelfiles <config.yaml>

Writes a finished fit's EXOFASTv2-compatible model files -- the model,
residuals and pretty model of every dataset, the SED residuals and
atmospheres and the MIST tracks, in EXOFASTv2's names, units and formats --
to ``<dirname(prefix)>/modelfiles/``, from the config it ran with and its
saved ``<prefix>_trace.nc``, without re-sampling.  Every fit already writes
them; this is for a fit that predates them, or for rewriting them.  The
files and their definitions: ``docs/model_files_plan.md``; the writer:
``outputs/model_files.create_model_files``.

Logging goes to the console only.  ``setup_logging`` would open
``<prefix>.log`` for writing, and that file is the fit's own log.
"""

import click

from .cli_summary import _console_logging
from .outputs.model_files import create_model_files


@click.command()
@click.argument("config_file")
@click.option(
    "--start",
    is_flag=True,
    help="Also write the start set, at the seeds the config writes.",
)
@click.option(
    "--logger-level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING"], case_sensitive=False),
    help="Console logging level.",
)
def main(config_file, start, logger_level):
    """Write the EXOFASTv2-compatible model files of a finished fit.

    CONFIG_FILE is the same system YAML passed to `exozippy`; its `prefix:`
    locates the saved trace (<prefix>_trace.nc), and the files go to the
    modelfiles/ directory beside it.
    """
    _console_logging(logger_level)
    written = create_model_files(config_file, start=start)
    if not written:
        click.echo("No component of this fit writes model files.")
        return
    click.echo(f"Wrote {len(written)} files to {written[-1].parent}")


if __name__ == "__main__":
    main()
