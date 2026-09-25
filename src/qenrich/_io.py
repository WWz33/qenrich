"""Object cache: parsed annotation "objects" stored as standard net TSVs."""

import gzip
import json
from pathlib import Path

import pandas as pd

from . import __version__


def open_text(path: str | Path):
    """Open a possibly-gzipped text file transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return open(path, encoding="utf-8-sig")


def read_names(path: str | Path) -> pd.DataFrame:
    """Read a multi-column id->names TSV (--desc) verbatim, dropping a header row.

    For ``go_zh.tsv`` (``ID\\tEnglish\\tChinese``) the caller picks column 2 as
    the English name and column 3 as the Chinese name. A first row whose first
    cell is not a term id (a header, e.g. ``id\\tname\\tname_zh``) is skipped.
    """
    df = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False, index_col=False)
    if len(df) and df.iloc[0, 0].strip().lower() in {"id", "term", "go", "gene"}:
        df = df.iloc[1:].reset_index(drop=True)
    return df


def cache_dir_for(annot_path: str | Path) -> Path:
    """Cache directory next to the annotation file: ``<file>.qenrich/``."""
    p = Path(annot_path)
    return p.with_name(p.name + ".qenrich")


def save_objects(objects: dict[str, pd.DataFrame], cdir: str | Path, source: str | Path, fmt: str) -> None:
    """Write each object as ``<cdir>/<name>.tsv`` plus a ``meta.json`` stamp."""
    cdir = Path(cdir)
    cdir.mkdir(parents=True, exist_ok=True)
    for name, df in objects.items():
        df.to_csv(cdir / f"{name}.tsv", sep="\t", index=False)
    meta = {
        "source": str(source),
        "mtime": Path(source).stat().st_mtime,
        "format": fmt,
        "version": __version__,
    }
    (cdir / "meta.json").write_text(json.dumps(meta))


def load_objects(cdir: str | Path) -> dict[str, pd.DataFrame]:
    """Read back every ``*.tsv`` object in a cache/db directory."""
    cdir = Path(cdir)
    objects = {p.stem: pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, index_col=False) for p in sorted(cdir.glob("*.tsv"))}
    if not objects:
        raise FileNotFoundError(f"no .tsv objects in {cdir}")
    return objects


def cache_fresh(cdir: str | Path, source: str | Path, fmt: str | None = None) -> bool:
    """True when the meta stamp matches the source mtime, the format and this qenrich.

    ``fmt`` is the format the caller is about to parse with: ``--format`` must not
    be defeated by a cache written from a different parse. A stamp written by an
    older qenrich (no ``version`` key) counts as stale, so upgrading the package
    re-parses instead of serving objects the current code would no longer emit.
    """
    meta = Path(cdir) / "meta.json"
    source = Path(source)
    if not (meta.is_file() and source.is_file()):
        return False
    try:
        stamp = json.loads(meta.read_text())
        fresh = stamp["mtime"] == source.stat().st_mtime and stamp["version"] == __version__
    except (KeyError, json.JSONDecodeError):
        return False
    return fresh and (fmt is None or stamp.get("format") == fmt)
