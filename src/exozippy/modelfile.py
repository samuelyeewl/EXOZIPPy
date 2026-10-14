"""
EXOFASTv2-compatible model files: the container a component returns.

EXOFASTv2 writes the model, residuals and "pretty" model of every dataset as
plain text in a ``modelfiles/`` subdirectory, so users can draw their own
figures without rebuilding the model.  ``Component.model_files(system,
point)`` returns a list of :class:`ModelFile`, one per file, and
``outputs/model_files.py`` writes them.

The split follows the architecture rule: WHICH quantity a file holds (the
full model, one planet's term, the trend), in WHICH unit and number format,
is component knowledge, so it lives in the component; placing the files --
directory, names, the index, stale-file removal -- is generic and lives in
the writer, which names no component.  The file definitions and why they
are what they are: ``docs/model_files_plan.md``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

#: EXOFASTv2's number formats, as IDL format codes translate: ``f0.N`` is
#: ``%.Nf`` (minimal width, N decimals), ``e0.N`` is ``%.Ne``, and ``aW`` is a
#: W-wide field, right-justified and truncated to W characters as Fortran's
#: (and IDL's) ``A`` edit descriptor does on output.
F6 = "%.6f"
F8 = "%.8f"
F10 = "%.10f"
F5 = "%.5f"
E6 = "%.6e"


def a_field(width):
    """The printf spelling of IDL's ``aW`` output format."""
    return f"%{width}.{width}s"


@dataclass
class ModelFile:
    """One EXOFASTv2-style text file.

    Parameters
    ----------
    kind : str
        The name between the fit's base and the extension, e.g.
        ``"model.transit_000"``; the file is ``<base>.<kind><ext>``.
    columns : sequence
        One 1-D array (or sequence of strings, for a text column) per
        column, all the same length.
    formats : sequence of str
        One printf code per column (``F8``, ``E6``, ``a_field(20)``, ...).
        Columns are joined by one space, as IDL's ``x`` separator writes.
    comment : str, optional
        A single header line, written verbatim.  ``None`` -- the default,
        and EXOFASTv2's ``/nocomment`` -- writes no header.
    ext : str
        ``".txt"`` except where EXOFASTv2 uses another (``".track"``).
    dataset : (str, str), optional
        ``(id, description)`` for the fit's index file: ``id`` the
        EXOFASTv2-style dataset name the file is numbered by
        (``"transit_000"``), ``description`` the dataset's name and data
        file.  EXOFASTv2's files carry no header, so the index is what says
        which file ``transit_000`` is.
    """

    kind: str
    columns: Sequence
    formats: Sequence[str]
    comment: Optional[str] = None
    ext: str = ".txt"
    dataset: Optional[tuple] = field(default=None)

    def lines(self):
        """The file's lines (no trailing newlines), header included.

        RAISES on a column count that differs from the format count, or on
        columns of unequal length: both are a bug in the component that
        built the file, never something to pad or truncate.
        """
        if len(self.columns) != len(self.formats):
            raise ValueError(
                f"model file {self.kind!r}: {len(self.columns)} columns but "
                f"{len(self.formats)} formats."
            )
        cols = [
            np.asarray(c) if not _is_text(c) else list(c) for c in self.columns
        ]
        lengths = {len(c) for c in cols}
        if len(lengths) > 1:
            raise ValueError(
                f"model file {self.kind!r}: columns of unequal length "
                f"{[len(c) for c in cols]}."
            )
        out = [] if self.comment is None else [self.comment]
        n = lengths.pop() if lengths else 0
        for r in range(n):
            out.append(
                " ".join(
                    fmt % _cell(c[r]) for fmt, c in zip(self.formats, cols)
                )
            )
        return out


def _is_text(column):
    return len(column) > 0 and isinstance(column[0], str)


def _cell(value):
    """A numpy scalar as the Python type printf formats natively."""
    if isinstance(value, str):
        return value
    return float(value)
