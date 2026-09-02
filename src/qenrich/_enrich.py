"""Bridge gene sets to decoupler's ORA (Fisher exact + BH FDR) and GSEA."""

import re

import numpy as np
import pandas as pd
import scipy.stats as sts

import decoupler as dc


def run_ora(
    net: pd.DataFrame,
    sets: dict[str, list[str]],
    tmin: int = 5,
    bg: list[str] | None = None,
    verbose: bool = False,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, dict[str, dict]]:
    """Enrich every gene set against a net via ``dc.mt.ora``.

    Each set becomes one row of a binary matrix over the universe (all net
    targets), ``n_up=len(set)`` recovers exactly the member genes, and the
    feature-specific background makes the Fisher table use the whole universe.

    Returns
    -------
    results : dict[str, pd.DataFrame]
        Per-set result table: term, term_size, overlap, genes, pvalue, log_or, padj.
    es_wide : pd.DataFrame
        set x term log odds ratios (for barplots).
    stats : dict[str, dict]
        Per-set run statistics (n_input, n_universe_hit, n_terms, n_pruned).
    """
    if bg is not None:
        bgset = set(bg)
        net = net[net["target"].isin(bgset)]
        if net.empty:
            raise ValueError("no overlap between annotation and --bg genes")
    universe = pd.Index(sorted(set(net["target"])))
    term_targets = net.groupby("source")["target"].apply(lambda s: set(s)).to_dict()
    results: dict[str, pd.DataFrame] = {}
    stats: dict[str, dict] = {}
    es_rows: dict[str, pd.Series] = {}
    for name, genes in sets.items():
        gs = sorted(set(genes) & set(universe))
        n_input = len(set(genes))
        if not gs:
            results[name] = pd.DataFrame(columns=["term", "term_size", "overlap", "genes", "pvalue", "log_or", "padj"])
            stats[name] = {"n_input": n_input, "n_hit": 0, "n_terms": 0, "n_pruned": len(term_targets)}
            continue
        row = pd.DataFrame(
            np.isin(universe, gs).astype(float).reshape(1, -1),
            index=[name],
            columns=universe,
        )
        es, pv = dc.mt.ora(row, net, n_up=len(gs), n_bm=0, n_bg=None, tmin=tmin, empty=False, verbose=verbose)
        terms = es.columns
        n, big_n = len(gs), len(universe)
        recs = []
        for t in terms:
            tg = term_targets[t]
            k = len(tg & set(gs))
            genes_hit = ";".join(sorted(tg & set(gs)))
            # two-sided Fisher raw p, consistent with decoupler's BH-adjusted padj
            tab = [[k, n - k], [len(tg) - k, big_n - len(tg) - n + k]]
            pval = sts.fisher_exact(tab, alternative="two-sided")[1]
            recs.append((t, len(tg), k, genes_hit, pval, float(es.loc[name, t]), float(pv.loc[name, t])))
        df = pd.DataFrame(recs, columns=["term", "term_size", "overlap", "genes", "pvalue", "log_or", "padj"])
        df = df.sort_values("padj").reset_index(drop=True)
        results[name] = df
        es_rows[name] = es.loc[name]
        stats[name] = {"n_input": n_input, "n_hit": len(gs), "n_terms": len(terms), "n_pruned": len(term_targets) - len(terms)}
    es_wide = pd.DataFrame(es_rows).T.reindex(columns=sorted({t for df in results.values() for t in df["term"]}))
    return results, es_wide, stats


def strip_suffix(sets: dict[str, list[str]], numeric: dict[str, dict[str, float]], net: pd.DataFrame):
    """Drop ``.1``-style version suffixes from genes and net targets (--strip-suffix)."""
    sets = {k: [re.sub(r"\.\d+$", "", g) for g in v] for k, v in sets.items()}
    numeric = {k: {re.sub(r"\.\d+$", "", g): w for g, w in v.items()} for k, v in numeric.items()}
    net = net.copy()
    net["target"] = net["target"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    net = net.drop_duplicates(subset=["source", "target"]).reset_index(drop=True)
    return sets, numeric, net


def _leading_edge(vec: dict[str, float], term_targets: set[str]) -> tuple[int, list[str]]:
    """Count and members of the leading-edge subset (clusterProfiler core_enrichment).

    Replicates decoupler's running sum: genes ranked by descending weight, hits
    add |w|/sum(|w in set|), misses subtract 1/(N-n); the running-sum peak cuts
    the leading edge.
    """
    genes = sorted(set(vec) & term_targets)
    order = sorted(vec, key=lambda g: -vec[g])
    n_all, n_set = len(order), len(genes)
    if not genes:
        return 0, []
    if n_all == n_set:  # term covers all ranked genes; all are leading edge
        return n_set, genes
    sum_set = sum(abs(vec[g]) for g in genes) or 1.0
    dec = 1.0 / (n_all - n_set)
    cum = 0.0
    mx_pos = mx_neg = 0.0
    j_pos = j_neg = -1
    for rank, g in enumerate(order):
        if g in genes:
            cum += abs(vec[g]) / sum_set
            if cum > mx_pos:
                mx_pos, j_pos = cum, rank
        else:
            cum -= dec
            if cum < mx_neg:
                mx_neg, j_neg = cum, rank
    if mx_pos > -mx_neg:
        edge = [g for g in order[: j_pos + 1] if g in genes]
    else:
        edge = [g for g in order[j_neg + 1:] if g in genes]
    return len(edge), edge


def run_gsea(
    net: pd.DataFrame,
    numeric: dict[str, dict[str, float]],
    tmin: int = 5,
    verbose: bool = False,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, dict[str, dict]]:
    """Enrich ranked vectors (gene -> weight, e.g. log2FC) via ``dc.mt.gsea``.

    Each vector becomes one row over the genes it contains; genes absent from
    the vector are not ranked (standard GSEA semantics).

    Returns per-set tables (term, term_size, nes, padj), a set x term NES wide
    frame, and run statistics.
    """
    term_targets = net.groupby("source")["target"].apply(lambda s: set(s)).to_dict()
    universe = set(net["target"])
    results: dict[str, pd.DataFrame] = {}
    stats: dict[str, dict] = {}
    es_rows: dict[str, pd.Series] = {}
    for name, vec in numeric.items():
        genes = sorted(set(vec) & universe)
        stats[name] = {"n_input": len(vec), "n_hit": len(genes), "n_terms": 0, "n_pruned": 0}
        if not genes:
            results[name] = pd.DataFrame(columns=["term", "term_size", "Count", "genes", "nes", "padj", "GeneRatio"])
            continue
        if all(w == 0 for w in (vec[g] for g in genes)):
            # all-zero weights crash decoupler's running sum (division by zero)
            results[name] = pd.DataFrame(columns=["term", "term_size", "Count", "genes", "nes", "padj", "GeneRatio"])
            continue
        if len(genes) < tmin:
            results[name] = pd.DataFrame(columns=["term", "term_size", "Count", "genes", "nes", "padj", "GeneRatio"])
            continue
        row = pd.DataFrame(
            [[vec[g] for g in genes]],
            index=[name],
            columns=genes,
        )
        try:
            es, pv = dc.mt.gsea(row, net, tmin=tmin, empty=False, verbose=verbose)
        except AssertionError as e:  # no term shares >= tmin targets with this row
            print(f"[qenrich] warning: GSEA set '{name}' skipped ({e})", file=sys.stderr)
            results[name] = pd.DataFrame(columns=["term", "term_size", "Count", "genes", "nes", "padj", "GeneRatio"])
            continue
        recs = []
        for t in es.columns:
            count, edge = _leading_edge({g: vec[g] for g in genes}, term_targets[t])
            recs.append((t, len(term_targets[t]), count, ";".join(sorted(edge)),
                         float(es.loc[name, t]), float(pv.loc[name, t])))
        df = pd.DataFrame(recs, columns=["term", "term_size", "Count", "genes", "nes", "padj"])
        df["padj"] = df["padj"].clip(lower=np.finfo(float).eps)  # 0.0 breaks log colour scales
        df["GeneRatio"] = df["Count"] / max(len(genes), 1)
        df = df.sort_values("padj").reset_index(drop=True)
        results[name] = df
        es_rows[name] = es.loc[name]
        stats[name]["n_terms"] = len(es.columns)
        stats[name]["n_pruned"] = len(term_targets) - len(es.columns)
    nes_wide = pd.DataFrame(es_rows).T.reindex(columns=sorted({t for df in results.values() for t in df["term"]}))
    return results, nes_wide, stats


def drop_parents(results: dict[str, pd.DataFrame], go, thr: float = 0.05) -> dict[str, pd.DataFrame]:
    """Collapse GO results: drop a parent term when a significant child exists."""
    out = {}
    for name, df in results.items():
        sig = set(df.loc[df["padj"] < thr, "term"])
        keep = [t for t in df["term"] if not (go.children(t) & sig)]
        out[name] = df[df["term"].isin(keep)].reset_index(drop=True)
    return out
