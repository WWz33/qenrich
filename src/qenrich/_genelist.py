"""Read gene list files where each column is one gene set (or a ranked vector)."""

import re
from pathlib import Path

from ._io import open_text

_SKIP = {"", "-", "NA", "N/A", "nan", "None"}
_NUM_CELL = re.compile(r"^(\S+?)[,;]([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)$")


def _split_line(line: str) -> list[str]:
    """Split on tabs if present, else on any whitespace."""
    if "\t" in line:
        return line.rstrip("\r\n").split("\t")
    return line.split()


def _has_header(rows: list[list[str]]) -> bool:
    """Treat the first row as a header unless any of its cells re-appears below.

    Gene identifiers almost always repeat across set columns (a gene belongs to
    several sets) or within the id column, while set names do not. Files where
    every gene is unique will be misdetected; ``--no-header`` forces the other
    way.
    """
    if len(rows) < 2:
        return False
    first = {c.strip() for c in rows[0] if c.strip() and c.strip() not in _SKIP}
    if not first:
        return False
    for row in rows[1:]:
        for c in row:
            if c.strip() and c.strip() in first:
                return False
    return True


def _is_numeric(cells: list[str]) -> bool:
    """A column is numeric (a ranked vector) when ~all cells are ``gene,weight``."""
    vals = [c.strip() for c in cells if c.strip() and c.strip() not in _SKIP]
    if not vals:
        return False
    hits = sum(1 for c in vals if _NUM_CELL.match(c))
    return hits > 0 and hits >= len(vals) - 1


def _parse_numeric(cells: list[str]) -> dict[str, float]:
    out = {}
    for c in cells:
        m = _NUM_CELL.match(c.strip())
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def read_genelist(
    path: str | Path, no_header: bool = False
) -> tuple[list[str], dict[str, list[str]], dict[str, dict[str, float]]]:
    """Parse a gene list file: one column per gene set.

    Parameters
    ----------
    path : str | Path
        Whitespace/tab separated file. Empty cells and ``-``/``NA`` are skipped.
    no_header : bool, optional
        Assume no header row and name sets ``set1..setN``. Default (False)
        auto-detects; falls back to ``set1..setN`` when no header is found.

    Returns
    -------
    columns : list[str]
        The detected column names.
    sets : dict[str, list[str]]
        Plain gene sets (for ORA), ordered and de-duplicated.
    numeric : dict[str, dict[str, float]]
        Columns whose cells look like ``gene,weight`` (for GSEA), gene -> weight.
    """
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            cells = _split_line(line)
            while cells and cells[-1] == "":
                cells.pop()  # trailing-tab artifact: empty last cell carries no data
            if cells:
                rows.append(cells)
    if not rows:
        raise ValueError(f"gene list file is empty: {path}")
    ncols = max(len(r) for r in rows)
    rows = [r + [""] * (ncols - len(r)) for r in rows]
    if no_header:
        header = [f"set{i + 1}" for i in range(ncols)]
        data = rows
    elif len(rows) == 1:
        # single-row file: almost always a header with no genes yet
        header = [c.strip() or f"set{i + 1}" for i, c in enumerate(rows[0])]
        data = []
    elif _has_header(rows):
        header = [c.strip() or f"set{i + 1}" for i, c in enumerate(rows[0])]
        data = rows[1:]
    else:
        header = [f"set{i + 1}" for i in range(ncols)]
        data = rows
    if len(set(header)) != len(header):
        raise ValueError(f"duplicate set names in gene list header: {header}")
    sets: dict[str, list[str]] = {}
    numeric: dict[str, dict[str, float]] = {}
    for j, name in enumerate(header):
        col = [r[j] for r in data]
        if _is_numeric(col):
            numeric[name] = _parse_numeric(col)
            continue
        seen: list[str] = []
        for c in col:
            g = c.strip()
            if g in _SKIP or g in seen:
                continue
            seen.append(g)
        sets[name] = seen
    return header, sets, numeric
