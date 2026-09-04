"""enrichplot (clusterProfiler) style plots, replicated from the R sources.

Visual spec extracted from enrichplot/R (dotplot.R, barplot.R, heatplot.R,
color_utils.R) and ggfun::theme_dose:

  - dotplot (.dotplot_internal): x = GeneRatio, y = Description (ordered by
    GeneRatio desc, labels truncated at 30 chars), shape-21 points (black
    border, fill = p.adjust), fill scale = two-colour gradient
    #327eba (low/significant) -> #e06663 (high) with log10 transform
    (get_enrichplot_color(2) reversed), size = Count via
    scale_size(range=c(3,8)) with pretty breaks, theme_dose.
  - barplot: x = Count, same fill scale and theme.
  - heatplot (no foldChange): black shape-21 dots on a white ground at
    gene x term membership, y = term labels, x = genes.
"""

from pathlib import Path

import numpy as np
import pandas as pd

SIG_LOW = "#327eba"   # smallest p.adjust
SIG_HIGH = "#e06663"  # largest p.adjust
LAB_MAX = 0           # 0 = no truncation (full term names)
FONT = 14             # theme_dose(font.size=14)


def _trunc(s: str) -> str:
    return s if LAB_MAX <= 0 or len(s) <= LAB_MAX else s[:LAB_MAX] + "..."


def _theme_dose(ax, title: str) -> None:
    """ggfun::theme_dose(font.size=14): theme_bw minus grid, bold centred title."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ax.set_facecolor("white")
    for side in ax.spines.values():
        side.set_visible(True)
        side.set_linewidth(0.8)
        side.set_color("black")
    ax.grid(False)
    ax.set_title(title, fontsize=FONT, fontweight="bold", loc="center", pad=10)
    ax.tick_params(axis="both", labelsize=int(FONT * 0.8), length=3)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_color("black")


def _sig_scale():
    """Fill colours + log10 normalisation replicating set_enrichplot_color defaults."""
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    cmap = LinearSegmentedColormap.from_list("enrichplot", [SIG_LOW, SIG_HIGH])
    return cmap, LogNorm


def _padj_norm(padj: np.ndarray):
    from matplotlib.colors import LogNorm

    # floor at 1e-16: permutation p-values can underflow to 0.0, which would
    # stretch a LogNorm to 1e-300 and pin every point to one end
    lo = max(float(np.nanmin(padj)), 1e-16)
    hi = max(float(np.nanmax(padj)), lo)
    if hi <= lo * 1.001:  # single unique padj: centre it like ggplot's scale expansion
        lo, hi = lo / 100, hi * 100
    return LogNorm(vmin=lo, vmax=hi)


def _size_map(counts: np.ndarray):
    """scale_size(range=c(3,8)): Count -> ggplot size 3..8 -> marker area pt^2."""
    lo, hi = float(np.nanmin(counts)), float(np.nanmax(counts))
    if hi == lo:
        hi = lo + 1
    def to_s(c):
        g = 3 + 5 * (np.asarray(c) - lo) / (hi - lo)
        return (g * 1.55) ** 2  # ggplot diameter -> matplotlib area
    breaks = np.unique(np.geomspace(max(lo, 1), max(hi, 2), 4).astype(int))
    return to_s, breaks


def _prep(df: pd.DataFrame, stats: dict, name_of, top: int) -> pd.DataFrame:
    d = df.copy()
    d["Description"] = [name_of(t) or t for t in d["term"]] if name_of else d["term"]
    d["Description"] = d["Description"].map(_trunc)
    if "Count" not in d.columns:  # ORA tables: Count = overlap
        d["Count"] = d["overlap"]
    if "GeneRatio" not in d.columns:
        d["setSize"] = max(stats.get("n_hit", 0), 1)
        d["GeneRatio"] = d["Count"] / d["setSize"]
    return d.sort_values("GeneRatio", ascending=False).head(top)


def _dotplot(d: pd.DataFrame, path: Path, title: str, figsize=(6.0, 4.0)) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    d = d.iloc[::-1]  # highest GeneRatio at top
    cmap, _ = _sig_scale()
    norm = _padj_norm(d["padj"].to_numpy())
    to_s, size_breaks = _size_map(d["Count"].to_numpy())
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    ax.scatter(
        d["GeneRatio"], d["Description"],
        s=to_s(d["Count"]),
        c=cmap(norm(d["padj"].to_numpy())),
        edgecolors="black", linewidths=0.6, zorder=3,
    )
    # Count size legend (scale_size pretty breaks)
    handles = [
        Line2D([], [], marker="o", linestyle="", markerfacecolor="none",
               markeredgecolor="black", markersize=np.sqrt(to_s(b)) / 1.0, label=str(b))
        for b in size_breaks
    ]
    leg = ax.legend(handles=handles, title="Count", loc="center left",
                    bbox_to_anchor=(1.02, 0.7), frameon=False,
                    labelspacing=1.4, borderpad=0.6, handletextpad=0.8)
    leg.get_title().set_fontsize(int(FONT * 0.85))
    for t in leg.get_texts():
        t.set_fontsize(int(FONT * 0.7))
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.7, pad=0.14,
    )
    cbar.set_label("p.adjust", size=int(FONT * 0.85))
    cbar.ax.tick_params(labelsize=int(FONT * 0.7))
    ax.set_xlabel("GeneRatio", fontsize=FONT)
    ax.set_ylabel(None)
    _theme_dose(ax, title)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _barplot(d: pd.DataFrame, path: Path, title: str, figsize=(6.0, 4.0)) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable

    d = d.sort_values("Count", ascending=True)  # highest Count at top
    cmap, _ = _sig_scale()
    norm = _padj_norm(d["padj"].to_numpy())
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    ax.barh(d["Description"], d["Count"], color=cmap(norm(d["padj"].to_numpy())), height=0.7)
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.7, pad=0.06)
    cbar.set_label("p.adjust", size=int(FONT * 0.85))
    cbar.ax.tick_params(labelsize=int(FONT * 0.7))
    ax.set_xlabel("Count", fontsize=FONT)
    ax.set_ylabel(None)
    _theme_dose(ax, title)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _heatplot(results_dfs: list[pd.DataFrame], path: Path, name_of, top_terms: int = 30, top_genes: int = 30) -> None:
    """No-foldChange heatplot: black shape-21 dots for gene x term membership."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows: list[tuple[str, str]] = []
    for d in results_dfs:
        for _, r in d.head(top_terms).iterrows():
            for g in [g for g in str(r["genes"]).split(";") if g]:
                rows.append((r["term"], g))
    if not rows:
        print(f"[qenrich] note: no gene memberships in the top terms, skipping {path.name}")
        return
    mat = pd.DataFrame(rows, columns=["term", "gene"]).drop_duplicates()
    terms = mat["term"].value_counts().head(top_terms).index
    mat = mat[mat["term"].isin(terms)]
    genes = mat["gene"].value_counts().head(top_genes).index
    mat = mat[mat["gene"].isin(genes)]
    piv = mat.assign(v=1).pivot(index="gene", columns="term", values="v")
    col_labels = [_trunc(name_of(t) or t) if name_of else _trunc(t) for t in piv.columns]
    fig, ax = plt.subplots(figsize=(1.0 + 0.16 * max(len(g) for g in piv.index) + 0.25 * len(col_labels),
                                    max(3.0, 0.26 * len(piv) + 1.0)), dpi=150)
    ax.scatter(
        np.repeat(np.arange(len(col_labels)), len(piv)),
        np.tile(np.arange(len(piv)), len(col_labels)),
        s=[60 * piv.iloc[i, j] for j in range(len(col_labels)) for i in range(len(piv))],
        color="black", marker="o",
    )
    ax.set_xticks(range(len(col_labels)), col_labels, rotation=45, ha="right", fontsize=int(FONT * 0.6))
    ax.set_yticks(range(len(piv)), piv.index, fontsize=int(FONT * 0.55))
    ax.set_xlim(-0.5, len(col_labels) - 0.5)
    ax.set_ylim(len(piv) - 0.5, -0.5)
    _theme_dose(ax, "")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_results_enrichplot(
    outdir: str | Path,
    results: dict[str, pd.DataFrame],
    stats: dict[str, dict],
    label_map: dict | None,
    top: int = 10,
    tag: str = "ep",
) -> None:
    """enrichplot-style per-set plots: GeneRatio dotplot, Count barplot, heatplot."""
    from ._plot import use_cjk_font

    use_cjk_font()
    name_of = (label_map or {}).get
    outdir = Path(outdir)
    prepped = {}
    for name, df in results.items():
        if not df.empty:
            prepped[name] = _prep(df, stats.get(name, {}), name_of, top)
    lab_len = max((d["Description"].str.len().max() for d in prepped.values() if len(d)), default=0)
    figsize = (max(5.5, 1.2 + 0.06 * lab_len), max(3.0, 0.28 * top + 1.0))
    for name, d in prepped.items():
        fname = name.replace("/", "_")
        _dotplot(d, outdir / f"{fname}_{tag}_dotplot.png", name, figsize)
        _barplot(d, outdir / f"{fname}_{tag}_barplot.png", name, figsize)
    # the tag keeps ORA and GSEA heatplots from overwriting each other
    _heatplot(list(results.values()), outdir / f"{tag}_heatplot.png", name_of, top_terms=top)
