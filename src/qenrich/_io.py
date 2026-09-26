"""Object cache: parsed annotation "objects" stored as standard net TSVs."""

import gzip
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__

# Bump when anything a stamp covers changes meaning (parse output, npz layout,
# stamp fields). Lets cache fixes invalidate old caches without a release bump.
CACHE_SCHEMA = 3

_PROP_STEM = "go_propagated"  # the propagated-GO cache, not an object

SKIP_CELLS = frozenset({"", "-", "NA", "N/A", "nan", "None"})


def atomic_write(path: str | Path, write, tmp_suffix: str = "") -> None:
    """Write via a temp sibling and ``os.replace``: a crash mid-write leaves the
    previous file intact instead of a half-written cache under the real name.

    ``tmp_suffix`` exists for writers that dispatch on the extension (np.savez
    only writes names ending in ".npz").
    """
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp" + tmp_suffix)
    write(tmp)
    os.replace(tmp, path)


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

    The stamp carries per-object term/gene counts (the run log reports them
    without reading objects a run does not use) and doubles as the ownership
    manifest: object TSVs it lists but this parse no longer produces are deleted,
    along with the propagated-GO cache built from the previous objects.
    """
    cdir = Path(cdir)
    cdir.mkdir(parents=True, exist_ok=True)
    meta_path = cdir / "meta.json"
    try:
        prev = json.loads(meta_path.read_text())
    except (OSError, json.JSONDecodeError):
        prev = {}
    for name in set(prev.get("objects", {})) - set(objects):
        try:
            (cdir / f"{name}.tsv").unlink()
        except OSError:
            pass
    for stale in (f"{_PROP_STEM}.npz", f"{_PROP_STEM}.json"):
        try:
            (cdir / stale).unlink()
        except OSError:
            pass
    for name, df in objects.items():
        atomic_write(cdir / f"{name}.tsv",
                     lambda p, df=df: df.to_csv(p, sep="\t", index=False))
    src = Path(source)
    meta = {
        "schema": CACHE_SCHEMA,
        "source": str(source),
        "mtime": src.stat().st_mtime,
        "size": src.stat().st_size,
        "format": fmt,
        "version": __version__,
        "objects": {name: {"terms": int(df["source"].nunique()),
                           "genes": int(df["target"].nunique())}
                    for name, df in objects.items()},
    }
    atomic_write(meta_path, lambda p: p.write_text(json.dumps(meta)))


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
    payload = dict(t=t_codes.astype("int32"), g=g_codes.astype("int32"),
                   tu=np.asarray(t_uniq, dtype=str), gu=np.asarray(g_uniq, dtype=str))
    atomic_write(path, lambda p: np.savez(p, **payload), tmp_suffix=".npz")


def load_net(path: str | Path) -> pd.DataFrame:
    """Read back what :func:`save_net` wrote, as ``string`` columns."""
    with np.load(path) as z:
        t, g, tu, gu = z["t"], z["g"], z["tu"], z["gu"]
    if len(t) != len(g):
        raise ValueError(f"inconsistent net cache {path}: {len(t)} vs {len(g)} rows")
    return pd.DataFrame({"source": pd.array(tu, dtype="string").take(t),
                         "target": pd.array(gu, dtype="string").take(g)})


def cache_fresh(cdir: str | Path, source: str | Path, fmt: str | None = None) -> bool:
    """True when the meta stamp matches the source (mtime and size), the format,
    the cache schema and this qenrich.

    ``fmt`` is the format the caller is about to parse with: ``--format`` must not
    be defeated by a cache written from a different parse. A stamp missing any
    checked field (older qenrich, older schema) counts as stale, so upgrading
    re-parses instead of serving objects the current code would no longer emit.
    """
    meta = Path(cdir) / "meta.json"
    source = Path(source)
    if not (meta.is_file() and source.is_file()):
        return False
    try:
        stamp = json.loads(meta.read_text())
        fresh = (stamp.get("schema") == CACHE_SCHEMA
                 and stamp["mtime"] == source.stat().st_mtime
                 and stamp.get("size") == source.stat().st_size
                 and stamp["version"] == __version__)
    except (KeyError, json.JSONDecodeError, OSError):
        return False
    return fresh and (fmt is None or stamp.get("format") == fmt)


def safe_names(names) -> dict[str, str]:
    """Map set names to file-safe basenames, unique across the list.

    ``/`` becomes ``_``; that mapping is not injective ("a/b" vs "a_b"), so a
    colliding name gets ``#2``, ``#3``, ... instead of overwriting the other
    set's files. Callers must pass the same name list for a run's TSVs and plots
    so both agree on one file per set.
    """
    out: dict[str, str] = {}
    used: set[str] = set()
    for name in names:
        base = str(name).replace("/", "_")
        cand, n = base, 1
        while cand in used:
            n += 1
            cand = f"{base}#{n}"
        used.add(cand)
        out[name] = cand
    return out
