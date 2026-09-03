"""Object cache: parsed annotation "objects" stored as standard net TSVs."""

import gzip
import json
from pathlib import Path

import pandas as pd


def open_text(path: str | Path):
    """Open a possibly-gzipped text file transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return open(path, encoding="utf-8-sig")


def read_names(path: str | Path) -> pd.DataFrame:
    """Read a multi-column id->names TSV (--desc) verbatim.

    For ``go_zh.tsv`` (``ID\\tEnglish\\tChinese``) the caller picks column 2 as
    the English name and column 3 as the Chinese name.
    """
    return pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False, index_col=False)


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
    (cdir / "meta.json").write_text(json.dumps({"source": str(source), "mtime": Path(source).stat().st_mtime, "format": fmt}))


def load_objects(cdir: str | Path) -> dict[str, pd.DataFrame]:
    """Read back every ``*.tsv`` object in a cache/db directory."""
    cdir = Path(cdir)
    objects = {p.stem: pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, index_col=False) for p in sorted(cdir.glob("*.tsv"))}
    if not objects:
        raise FileNotFoundError(f"no .tsv objects in {cdir}")
    return objects


def cache_fresh(cdir: str | Path, source: str | Path) -> bool:
    """True when meta stamp matches current source file mtime."""
    meta = Path(cdir) / "meta.json"
    source = Path(source)
    if not (meta.is_file() and source.is_file()):
        return False
    try:
        return json.loads(meta.read_text())["mtime"] == source.stat().st_mtime
    except (KeyError, json.JSONDecodeError):
        return False
