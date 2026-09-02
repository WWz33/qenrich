"""Minimal go-basic.obo parser: term metadata, alt_id mapping, ancestor propagation."""

import re
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
        with open_text(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line == "[Term]":
                    cur = {"id": None, "name": None, "namespace": None, "is_a": set(), "part_of": set(), "alt": [], "obsolete": False}
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
                elif line == "" and cur is not None and cur["id"]:
                    if not cur["obsolete"]:
                        parents[cur["id"]] = cur["is_a"] | cur["part_of"]
                        meta[cur["id"]] = (cur["name"] or "", cur["namespace"] or "")
                        for a in cur["alt"]:
                            alt[a] = cur["id"]
                    cur = None
        if cur is not None and cur["id"] and not cur["obsolete"]:
            parents[cur["id"]] = cur["is_a"] | cur["part_of"]
            meta[cur["id"]] = (cur["name"] or "", cur["namespace"] or "")
            for a in cur["alt"]:
                alt[a] = cur["id"]
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

        alt_id entries are translated to their primary ID first.
        """
        rows = []
        for term, gene in zip(net["source"], net["target"], strict=True):
            term = self._alt.get(term, term)  # translate alt_id -> primary
            if term not in self._meta:
                continue  # drop unknown/obsolete terms
            rows.append((term, gene))
            for anc in self.ancestors(term):
                rows.append((anc, gene))
        out = pd.DataFrame(rows, columns=["source", "target"]).drop_duplicates()
        return out[out["source"].isin(self._meta)].reset_index(drop=True)
