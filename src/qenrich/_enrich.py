"""ORA over a term/gene net: a hypergeometric test per term, BH across the tested terms."""

import re

import numpy as np
import pandas as pd

from ._fisher import fisher_pvalues


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (the standard step-up, over the tested terms)."""
    p = np.asarray(p, dtype=float)
    m = p.size
    if m == 0:
        return p
    order = np.argsort(p, kind="stable")
    ranked = p[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(ranked, 0, 1)
    return out


def run_ora(
    net: pd.DataFrame,
    sets: dict[str, list[str]],
    tmin: int = 5,
    bg: list[str] | None = None,
    alternative: str = "greater",
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, dict[str, dict]]:
    """Enrich every gene set against the net (hypergeometric test + BH FDR).

    Reproduces clusterProfiler's ``enrichGO``: the p-value is the one-sided
    over-representation test ``phyper(k - 1, M, N - M, n, lower.tail = FALSE)``,
    computed by :mod:`qenrich._fisher`, and BH runs over the terms that hold at
    least one query gene, which are also the terms reported. Gene sets with no
    query gene are counted in the statistics but not tested.

    ``alternative='less'`` tests depletion instead (P(X <= k)); there is no
    clusterProfiler counterpart for it.

    Returns
    -------
    results : dict[str, pd.DataFrame]
        Per-set result table: term, term_size, overlap, genes, pvalue, log_or, padj.
    es_wide : pd.DataFrame
        set x term log odds ratios (for barplots).
    stats : dict[str, dict]
        Per-set run statistics (n_input, n_universe_hit, n_terms, n_pruned, n_empty).
    """
    if bg is not None:
        bgset = set(bg)
        net = net[net["target"].isin(bgset)]
        if net.empty:
            raise ValueError("no overlap between annotation and --bg genes")
    universe = pd.Index(sorted(pd.unique(net["target"])))
    big_n = len(universe)
    gene_pos = {g: i for i, g in enumerate(universe)}
    universe_set = set(universe)
    # One integer code per (term, gene) pair, deduplicated and sorted. Sorted codes
    # turn each term's gene positions into a contiguous slice, so no per-term
    # Python set work is needed. (Iterating pandas' Arrow-backed string columns in
    # Python costs seconds at organism scale, hence pd.unique/pd.factorize and no
    # set(net["target"]).)
    gene_codes, gene_uniques = pd.factorize(net["target"])
    term_codes_raw, term_names = pd.factorize(net["source"])
    pair_pos_all = np.array([gene_pos[g] for g in gene_uniques], dtype=np.int64)[gene_codes]
    pairs = term_codes_raw.astype(np.int64) * big_n + pair_pos_all
    if pairs.size:  # an empty net has no pairs, but the loop below still runs
        pairs.sort()
        keep = np.empty(len(pairs), dtype=bool)
        keep[0] = True
        np.not_equal(pairs[1:], pairs[:-1], out=keep[1:])
        pairs = pairs[keep]
    pair_terms, pair_pos = np.divmod(pairs, big_n)
    term_sizes = np.bincount(pair_terms, minlength=len(term_names))
    # +1 so every term's slice is [starts[c], starts[c + 1]); the last boundary is
    # the pair count
    term_starts = np.searchsorted(pair_terms, np.arange(len(term_names) + 1))
    sized_codes = np.where(term_sizes >= tmin)[0]  # terms passing --tmin
    gene_names = universe.to_numpy()

    results: dict[str, pd.DataFrame] = {}
    stats: dict[str, dict] = {}
    es_rows: dict[str, pd.Series] = {}
    for name, genes in sets.items():
        gs_set = set(genes) & universe_set
        n_input = len(set(genes))
        n_pruned = len(term_names) - len(sized_codes)
        if not gs_set:
            results[name] = pd.DataFrame(columns=["term", "term_size", "overlap", "genes", "pvalue", "log_or", "padj"])
            stats[name] = {"n_input": n_input, "n_hit": 0, "n_terms": 0,
                           "n_pruned": n_pruned, "n_empty": len(sized_codes)}
            continue
        n = len(gs_set)
        mask = np.zeros(big_n, dtype=bool)
        mask[[gene_pos[g] for g in gs_set]] = True
        overlap_all = np.bincount(pair_terms, weights=mask[pair_pos], minlength=len(term_names)).astype(np.int64)
        # clusterProfiler tests and reports only the gene sets holding a query gene
        tested_codes = sized_codes[overlap_all[sized_codes] > 0]
        tested_codes = tested_codes[np.argsort(term_names[tested_codes])]  # sorted(tested)
        sizes = term_sizes[tested_codes]
        observed = overlap_all[tested_codes]
        pvals = fisher_pvalues(
            sizes,
            np.full(len(tested_codes), n, dtype=np.int64),
            np.full(len(tested_codes), big_n, dtype=np.int64),
            observed,
            alternative=alternative,
        )
        recs = []
        for pos, code in enumerate(tested_codes):
            t, K, k = term_names[code], int(sizes[pos]), int(observed[pos])
            # the 2x2 table: a=k (both), b=term-only, c=set-only, d=neither
            a, b = k, K - k
            c, d = n - k, big_n - K - n + k
            lor = np.log((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)))  # Haldane-Anscombe
            in_term = pair_pos[term_starts[code]:term_starts[code + 1]]  # this term's genes
            genes_hit = ";".join(sorted(gene_names[in_term[mask[in_term]]]))
            recs.append((t, K, k, genes_hit, float(pvals[pos]), float(lor)))
        df = pd.DataFrame(recs, columns=["term", "term_size", "overlap", "genes", "pvalue", "log_or"])
        df["padj"] = _bh(df["pvalue"].values)
        # stable: integer tables repeat, so many terms tie exactly (same size and
        # overlap -> same p). Keeping the term order for ties makes the output
        # deterministic instead of whatever quicksort happens to produce.
        df = df.sort_values("padj", kind="stable").reset_index(drop=True)
        results[name] = df
        es_rows[name] = df.set_index("term")["log_or"].reindex(sorted(term_names[tested_codes])) if len(df) else pd.Series(dtype=float)
        stats[name] = {"n_input": n_input, "n_hit": len(gs_set), "n_terms": len(tested_codes),
                       "n_pruned": n_pruned, "n_empty": len(sized_codes) - len(tested_codes)}
    es_wide = pd.DataFrame(es_rows).T.reindex(columns=sorted({t for df in results.values() for t in df["term"]}))
    return results, es_wide, stats


def strip_suffix(sets: dict[str, list[str]], net: pd.DataFrame):
    """Drop ``.1``-style version suffixes from genes and net targets (--strip-suffix)."""
    sets = {k: [re.sub(r"\.\d+$", "", g) for g in v] for k, v in sets.items()}
    net = net.copy()
    net["target"] = net["target"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    net = net.drop_duplicates(subset=["source", "target"]).reset_index(drop=True)
    return sets, net


def drop_parents(results: dict[str, pd.DataFrame], go, thr: float = 0.05) -> dict[str, pd.DataFrame]:
    """Collapse GO results: drop a parent term when a significant child exists."""
    out = {}
    for name, df in results.items():
        sig = set(df.loc[df["padj"] < thr, "term"])
        keep = [t for t in df["term"] if not (go.children(t) & sig)]
        out[name] = df[df["term"].isin(keep)].reset_index(drop=True)
    return out
