"""Object cache: parsed annotation "objects" stored as standard net TSVs."""

import gzip
import json
from pathlib import Path

import pandas as pd

import numpy as np

from . import __version__

_PROP_STEM = "go_propagated"  # the propagated-GO cache, not an object

SKIP_CELLS = frozenset({"", "-", "NA", "N/A", "nan", "None"})


def open_text(path: str | Path):
    """Open a possibly-gzipped text file transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return open(path, encoding="utf-8-sig")


def read_names(path: str | Path) -> pd.DataFrame:
    """Read a multi-column id->names TSV (--zh TABLE) verbatim, dropping a header row.

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
    """Write each object as ``<cdir>/<name>.tsv`` plus a ``meta.json`` stamp.

    The stamp also carries per-object term/gene counts: the CLI reports them
    without reading objects a run does not use.
    """
    cdir = Path(cdir)
    cdir.mkdir(parents=True, exist_ok=True)
    for name, df in objects.items():
        df.to_csv(cdir / f"{name}.tsv", sep="\t", index=False)
    meta = {
        "source": str(source),
        "mtime": Path(source).stat().st_mtime,
        "format": fmt,
        "version": __version__,
        "objects": {name: {"terms": int(df["source"].nunique()),
                           "genes": int(df["target"].nunique())}
                    for name, df in objects.items()},
    }
    (cdir / "meta.json").write_text(json.dumps(meta))


def load_object(cdir: str | Path, name: str) -> pd.DataFrame:
    """Read one object TSV from a cache/db directory."""
    return pd.read_csv(Path(cdir) / f"{name}.tsv", sep="\t", dtype=str,
                       keep_default_na=False, index_col=False)


def load_objects(cdir: str | Path) -> dict[str, pd.DataFrame]:
    """Read back every ``*.tsv`` object in a cache/db directory.

    The propagated-GO cache is not an object: it is keyed to the raw ``go``
    object it was built from and read only by the CLI.
    """
    cdir = Path(cdir)
    objects = {p.stem: load_object(cdir, p.stem)
               for p in sorted(cdir.glob("*.tsv")) if p.stem != _PROP_STEM}
    if not objects:
        raise FileNotFoundError(f"no .tsv objects in {cdir}")
    return objects


def save_net(df: pd.DataFrame, path: str | Path) -> None:
    """Write a net as integer codes plus the distinct ids they index (npz).

    Reading a 1.8M-row net TSV back costs ~0.5s; the codes load in ~0.05s and
    rebuild the same string columns through Arrow ``take``.
    """
    t_codes, t_uniq = pd.factorize(df["source"])
    g_codes, g_uniq = pd.factorize(df["target"])
    np.savez(path, t=t_codes.astype("int32"), g=g_codes.astype("int32"),
             tu=np.asarray(t_uniq, dtype=str), gu=np.asarray(g_uniq, dtype=str))


def load_net(path: str | Path) -> pd.DataFrame:
    """Read back what :func:`save_net` wrote, as ``string`` columns."""
    z = np.load(path)
    return pd.DataFrame({"source": pd.array(z["tu"], dtype="string").take(z["t"]),
                         "target": pd.array(z["gu"], dtype="string").take(z["g"])})


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
