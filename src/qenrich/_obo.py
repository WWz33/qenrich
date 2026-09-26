"""Minimal go-basic.obo parser: term metadata, alt_id mapping, ancestor propagation."""

import pickle
import sys
from functools import lru_cache
from pathlib import Path

from ._io import CACHE_SCHEMA, atomic_write, open_text

import pandas as pd

from . import __version__


class GeneOntology:
    """Terms from an OBO file with true-path-rule propagation.

    Propagates along ``is_a`` and ``part_of`` only, the two relations that make
    up the go-basic backbone.
    """

    def __init__(self, parents: dict[str, set[str]], meta: dict[str, tuple[str, str]], alt: dict[str, str]):
        self._parents = parents
        self._meta = meta
        self._alt = alt
        self._children: dict[str, set[str]] = {}
        for child, ps in parents.items():
            for p in ps:
                self._children.setdefault(p, set()).add(child)

    @classmethod
    def from_obo(cls, path: str) -> "GeneOntology":
        parents: dict[str, set[str]] = {}
        meta: dict[str, tuple[str, str]] = {}
        alt: dict[str, str] = {}
        cur: dict[str, object] | None = None

        def commit(c: dict) -> None:
            """Register a stanza; obsolete terms only via their replacement."""
            if c["obsolete"]:
                # annotations to a retired term keep their signal: GO points at the
                # term that took over (replaced_by); consider= has no single target
                rep = c.get("replaced_by")
                if rep:
                    alt[c["id"]] = rep
                    for a in c["alt"]:
                        alt[a] = rep
                return
            parents[c["id"]] = c["is_a"] | c["part_of"]
            meta[c["id"]] = (c["name"] or "", c["namespace"] or "")
            for a in c["alt"]:
                alt[a] = c["id"]

        with open_text(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line == "[Term]":
                    cur = {"id": None, "name": None, "namespace": None, "is_a": set(),
                           "part_of": set(), "alt": [], "obsolete": False, "replaced_by": None}
                elif line.startswith("[") and line != "[Term]":
                    cur = None
                elif cur is not None and ": " in line:
                    key, val = line.split(": ", 1)
                    val = val.split("!")[0].strip()
                    if key == "id":
                        cur["id"] = val
                    elif key == "name":
                        cur["name"] = val
                    elif key == "namespace":
                        cur["namespace"] = val
                    elif key == "is_a":
                        cur["is_a"].add(val)
                    elif key == "relationship" and val.startswith("part_of "):
                        cur["part_of"].add(val.split()[1])
                    elif key == "alt_id":
                        cur["alt"].append(val)
                    elif key == "is_obsolete" and val == "true":
                        cur["obsolete"] = True
                    elif key == "replaced_by":
                        cur["replaced_by"] = val
                elif line == "" and cur is not None and cur["id"]:
                    commit(cur)
                    cur = None
        if cur is not None and cur["id"]:
            commit(cur)
        # a redirect must land on a live term: follow retired -> retired chains
        for k in list(alt):
            tgt, seen = alt[k], set()
            while tgt not in meta and tgt in alt and tgt not in seen:
                seen.add(tgt)
                tgt = alt[tgt]
            alt[k] = tgt
        return cls(parents, meta, alt)

    @classmethod
    def cached(cls, path: str, cache_dir: str | Path | None = None) -> "GeneOntology":
        """``from_obo``, reusing a pickle of the parsed stanzas when it is fresh.

        Parsing the 32MB go-basic.obo costs ~0.45s per run; the pickle is 4MB and
        loads in ~0.05s. The stamp records the cache schema, the OBO path, its
        mtime and this qenrich version, so a new OBO release, a new qenrich or a
        cache-format change re-parses. Without a cache dir (net TSV and --db
        inputs) there is nowhere to put it.
        """
        if cache_dir is None:
            return cls.from_obo(path)
        p = Path(cache_dir) / "obo.pkl"
        stamp = {"schema": CACHE_SCHEMA, "obo": str(Path(path).resolve()),
                 "obo_mtime": Path(path).stat().st_mtime, "version": __version__}
        try:
            with open(p, "rb") as fh:
                saved, payload = pickle.load(fh)
            if saved == stamp:
                return cls(*payload)
        except (OSError, EOFError, ValueError, TypeError, AttributeError, pickle.UnpicklingError):
            pass  # unreadable or stale pickle: parse and rewrite
        go = cls.from_obo(path)

        def dump(q):
            with open(q, "wb") as fh:
                pickle.dump((stamp, (go._parents, go._meta, go._alt)), fh, protocol=5)

        try:
            atomic_write(p, dump)
        except OSError as e:  # read-only cache dir must not abort the run
            print(f"[qenrich] warning: could not write {p} ({e}); continuing", file=sys.stderr)
        return go

    @lru_cache(maxsize=None)
    def ancestors(self, term: str) -> frozenset[str]:
        seen: set[str] = set()
        stack = list(self._parents.get(term, ()))
        while stack:
            t = stack.pop()
            if t in seen:
                continue
            seen.add(t)
            stack.extend(self._parents.get(t, ()))
        return frozenset(seen)

    def name(self, term: str) -> str:
        return self._meta.get(term, ("", ""))[0]

    def children(self, term: str) -> set[str]:
        return self._children.get(term, set())

    def namespace(self, term: str) -> str:
        return self._meta.get(term, ("", ""))[1]

    def propagate(self, net: pd.DataFrame) -> pd.DataFrame:
        """Expand a go net so every (gene, term) row repeats for all ancestors.

        alt_id entries are translated to their primary ID first. Terms with no
        entry in this ontology (e.g. an OBO older than the annotation file) are
        dropped; a warning reports the loss instead of hiding it.

        Vectorized: factorize puts the work on the ~40k distinct ids instead of
        the 660k-row Arrow string column (Python-level iteration over Arrow
        strings costs ~1s per column at this scale); genes are repeated with
        ``np.repeat`` and pairs are deduplicated through the integer codes.
        """
        import numpy as np

        raw_codes, raw_uniq = pd.factorize(net["source"])
        gene_codes, gene_uniq = pd.factorize(net["target"])
        alt, meta = self._alt, self._meta
        mapped = np.array([alt.get(t, t) for t in raw_uniq], dtype=object)
        alive = np.array([t in meta for t in mapped], dtype=bool)
        if not alive.any():
            return pd.DataFrame(columns=["source", "target"])
        anc = {t: sorted(self.ancestors(t)) for t in mapped[alive]}  # resolved once
        # per distinct source id: [term, *ancestors]; a dropped id expands to []
        per_term = [[t, *anc[t]] if t in anc else [] for t in mapped]
        # every row repeats its gene len(per_term) times: the term, then ancestors
        counts = np.fromiter((len(per_term[c]) for c in raw_codes), np.int64, len(raw_codes))
        # rows with a dropped term contribute nothing; every other row repeats
        # its own term's expansion [term, *ancestors]
        terms = np.concatenate([per_term[c] for c in raw_codes if per_term[c]])
        genes = np.repeat(gene_codes, counts)
        # dedupe (term, gene) pairs through integer codes
        tcodes, tuniq = pd.factorize(terms)
        pairs = tcodes.astype(np.int64) * len(gene_uniq) + genes
        pairs.sort()
        keep = np.empty(len(pairs), dtype=bool)
        keep[0] = True
        np.not_equal(pairs[1:], pairs[:-1], out=keep[1:])
        out = pd.DataFrame({"source": tuniq[pairs[keep] // len(gene_uniq)],
                            "target": gene_uniq[pairs[keep] % len(gene_uniq)]})
        out = out[out["source"].isin(meta)].reset_index(drop=True)
        dropped_rows = int((~alive[raw_codes]).sum())
        if dropped_rows:
            dropped_terms = {str(raw_uniq[i]) for i in np.where(~alive)[0]}
            # genes whose every annotation was dropped: compare unique-gene codes
            # against the survivors' codes (a Python string-set diff costs ~1s here)
            survived = pairs[keep] % len(gene_uniq)
            lost_mask = np.ones(len(gene_uniq), dtype=bool)
            lost_mask[survived] = False
            n_lost = int(lost_mask.sum())
            sample = ", ".join(sorted(dropped_terms)[:5])
            more = f" (+{len(dropped_terms) - 5} more)" if len(dropped_terms) > 5 else ""
            msg = (f"[qenrich] warning: {dropped_rows} annotation(s) dropped: "
                   f"{len(dropped_terms)} GO term(s) have no entry in the OBO "
                   f"(retired without a replacement, or the OBO predates the annotation): "
                   f"{sample}{more}")
            if n_lost:
                msg += f"; {n_lost} gene(s) left with no GO annotation"
            print(msg, file=sys.stderr)
        return out

