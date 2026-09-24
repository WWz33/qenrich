"""Minimal go-basic.obo parser: term metadata, alt_id mapping, ancestor propagation."""

import re
import sys
from functools import lru_cache

from ._io import open_text

import pandas as pd


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
        """
        rows = []
        dropped_terms: set[str] = set()
        dropped_rows = 0
        for raw, gene in zip(net["source"], net["target"], strict=True):
            term = self._alt.get(raw, raw)  # alt_id -> primary, retired -> replacement
            if term not in self._meta:
                dropped_terms.add(raw)  # the id as it appears in the user's file
                dropped_rows += 1
                continue
            rows.append((term, gene))
            for anc in self.ancestors(term):
                rows.append((anc, gene))
        out = pd.DataFrame(rows, columns=["source", "target"]).drop_duplicates()
        out = out[out["source"].isin(self._meta)].reset_index(drop=True)
        if dropped_terms:
            lost = set(net["target"]) - set(out["target"])
            sample = ", ".join(sorted(dropped_terms)[:5])
            more = f" (+{len(dropped_terms) - 5} more)" if len(dropped_terms) > 5 else ""
            msg = (f"[qenrich] warning: {dropped_rows} annotation(s) dropped: "
                   f"{len(dropped_terms)} GO term(s) have no entry in the OBO "
                   f"(retired without a replacement, or the OBO predates the annotation): "
                   f"{sample}{more}")
            if lost:
                msg += f"; {len(lost)} gene(s) left with no GO annotation"
            print(msg, file=sys.stderr)
        return out
