"""Read gene list files where each column is one gene set."""

from pathlib import Path

from ._io import SKIP_CELLS, open_text


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
    first = {c.strip() for c in rows[0] if c.strip() and c.strip() not in SKIP_CELLS}
    if not first:
        return False
    for row in rows[1:]:
        for c in row:
            if c.strip() and c.strip() in first:
                return False
    return True


def read_genelist(
    path: str | Path, no_header: bool = False
) -> tuple[list[str], dict[str, list[str]]]:
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
        Gene sets (for ORA), ordered and de-duplicated.
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
    for j, name in enumerate(header):
        seen: list[str] = []
        for r in data:
            c = r[j].strip()
            if c in SKIP_CELLS or c in seen:
                continue
            seen.append(c)
        sets[name] = seen
    return header, sets
