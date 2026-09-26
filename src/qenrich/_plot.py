"""Plotting helpers: per-set barplot/dotplot plus the summary heatmap.

The barplot and dotplot reproduce decoupler's ``dc.pl.barplot`` / ``dc.pl.dotplot``
(Copyright (c) 2025 scverse, Pau Badia i Mompel and Saez lab, BSD 3-Clause), with
the seaborn bar call and the ``Plotter`` wrapper rewritten in plain matplotlib so
that qenrich keeps only pandas/numpy/scipy/matplotlib as dependencies.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd


def _set_limits(vmin, vcenter, vmax, values: np.ndarray) -> tuple[float, float, float]:
    """Color limits for a diverging scale: mirror the range around ``vcenter``."""
    if vmin is None:
        vmin = values.min()
    if vmax is None:
        vmax = values.max()
    if vcenter is None:
        vcenter = values.mean()
    if vmin >= vcenter:
        vmin = -vmax
    if vcenter >= vmax:
        vmax = -vmin
    if vmin == vmax:  # every plotted score identical (e.g. all zero): keep the norm valid
        vmin, vmax = vmin - 1.0, vmax + 1.0
    return vmin, vcenter, vmax


def _unique_label_map(columns, label_map: dict) -> dict:
    """term -> label, guaranteed unique across ``columns``.

    Distinct terms can share a display name (GO has duplicate names); renaming
    them to the same label makes the plot average them into one bar/dot. On a
    collision the term id is used instead, and if that too is taken a numeric
    suffix is appended, so the returned labels are always distinct.
    """
    seen: set[str] = set()
    out: dict[str, str] = {}
    for c in columns:
        lab = label_map.get(c, c)
        if lab in seen:
            base = lab = c
            n = 2
            while lab in seen:
                lab = f"{base}#{n}"
                n += 1
        seen.add(lab)
        out[c] = lab
    return out


def _disp_len(s: str) -> float:
    """Display width in latin-char equivalents: CJK glyphs are ~1.7x wider."""
    return sum(1.7 if ord(c) > 0x2E80 else 1.0 for c in s)


def _barplot(fig, ax, data: pd.DataFrame, name: str, top: int = 15, cmap: str = "RdBu_r") -> None:
    """Horizontal bars of one set's top scores, faces colored by score.

    Top ``top`` terms by absolute score, drawn in descending order; the bar faces
    run through a zero-centered TwoSlopeNorm so enrichment and depletion read as
    different colors.
    """
    import matplotlib
    import matplotlib.pyplot as plt

    scores = data.loc[name].to_numpy(dtype=float)
    order = np.argsort(-np.abs(scores))[:top]          # top terms by |score|
    order = order[np.argsort(scores[order])][::-1]     # descending: largest first
    labels = data.columns.to_numpy()[order]
    values = scores[order]

    vmin, vcenter, vmax = _set_limits(None, 0, None, values)
    norm = matplotlib.colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    ax.barh(np.arange(len(values)), values, height=0.8, color=plt.get_cmap(cmap)(norm(values)))
    ax.set_yticks(np.arange(len(values)), labels)
    # seaborn's categorical axis: exact limits, first row at the top
    ax.set_ylim(-0.5, len(values) - 0.5)
    ax.invert_yaxis()
    ax.set_xlabel("Score")
    ax.set_ylabel("")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=0.5)


def _dotplot(fig, ax, df: pd.DataFrame, x: str, y: str, c: str, s: str,
             top: int = 15, scale: float = 0.15, cmap: str = "RdBu_r", vcenter=None) -> None:
    """Dots at ``x`` (score), sized by ``s`` (a -log10 padj scale), colored by ``c``.

    Top ``top`` rows by |x| drawn in ascending x, marker areas
    ``(s * scale * markersize)**2``, with a size legend and a colorbar titled by
    the column names.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    x_all = df[x].to_numpy(dtype=float)
    keep = np.argsort(-np.abs(x_all))[:top]
    keep = keep[np.argsort(x_all[keep])]
    x_vals, y_vals = x_all[keep], df[y].to_numpy()[keep]
    c_vals, s_vals = df[c].to_numpy(dtype=float)[keep], df[s].to_numpy(dtype=float)[keep]

    ns = (s_vals * scale * plt.rcParams["lines.markersize"]) ** 2
    ax.grid(axis="x")
    norm = TwoSlopeNorm(vmin=None, vcenter=vcenter, vmax=None) if vcenter is not None else None
    scatter = ax.scatter(x=x_vals, y=y_vals, c=c_vals, s=ns, cmap=cmap, norm=norm)
    ax.set_axisbelow(True)
    ax.set_xlabel(x)
    handles, labels = scatter.legend_elements(
        prop="sizes", num=3, fmt="{x:.2f}",
        func=lambda v: np.sqrt(v) / plt.rcParams["lines.markersize"] / scale,
    )
    ax.legend(handles, labels, title=s, frameon=False, loc="lower left", bbox_to_anchor=(1.05, 0.5),
              alignment="left", labelspacing=1.0)
    clb = fig.colorbar(scatter, ax=ax, shrink=0.25, aspect=5, orientation="vertical", anchor=(0.0, 0.0))
    clb.ax.set_title(c, loc="left")
    ax.margins(x=0.25, y=0.1)


def plot_results(outdir: str | Path, results: dict, es_wide: pd.DataFrame, label_map: dict | None) -> None:
    """Per set: a barplot of the scores and a dotplot of score against padj."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    use_cjk_font()  # term labels may be Chinese (--zh table.tsv)

    from ._io import safe_names

    outdir = Path(outdir)
    safe = safe_names(results)
    labels = _unique_label_map(es_wide.columns, label_map) if label_map else None
    if labels:
        es_wide = es_wide.rename(columns=labels)
    lab_len = max((_disp_len(str(c)) for c in es_wide.columns), default=0)
    figsize = (max(6.5, 3.0 + 0.11 * lab_len), 3.5)

    for name, df in results.items():
        if df.empty:
            continue
        fname = safe[name]
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=150, tight_layout=True)
        _barplot(fig, ax, es_wide.fillna(0.0), name)
        fig.savefig(outdir / f"{fname}_barplot.png", bbox_inches="tight")
        plt.close(fig)

        long = df.head(15).copy()
        # floor keeps non-significant terms visible
        long["nlp"] = (-np.log10(long["padj"].clip(lower=2.22e-16))).clip(lower=0.3)
        long["label"] = long["term"].map(labels).fillna(long["term"]) if labels else long["term"]
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=150, tight_layout=True)
        _dotplot(fig, ax, long, x="log_or", y="label", c="log_or", s="nlp")
        fig.savefig(outdir / f"{fname}_dotplot.png", bbox_inches="tight")
        plt.close(fig)


def plot_heatmap(summary: pd.DataFrame, outdir: str | Path, name_of) -> None:
    """term x set matrix of -log10(padj), top 30 terms."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    use_cjk_font()

    outdir = Path(outdir)
    p = summary.copy()
    p["nlp"] = -np.log10(p["padj"].clip(lower=2.22e-16))
    mat = p.pivot_table(index="term", columns="set", values="nlp", aggfunc="max")
    mat = mat.loc[mat.max(axis=1).sort_values(ascending=False).index].head(30)
    # no layout="tight" here: with ~30 rows of labels tight layout gives up and warns
    fig, ax = plt.subplots(figsize=(1.0 + 1.4 * mat.shape[1], 0.22 * len(mat) + 1.0), dpi=150)
    im = ax.imshow(mat.values, aspect="auto", cmap="Reds")
    ax.set_xticks(range(mat.shape[1]), mat.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(mat)), [name_of(t) or t for t in mat.index], fontsize=6)
    fig.colorbar(im, ax=ax, label="-log10(padj)", shrink=0.6)
    fig.savefig(outdir / "summary_heatmap.png", bbox_inches="tight")


_PROBE = "细胞过程Aa1"  # must cover CJK + latin + digits in ONE font


@lru_cache(maxsize=1)
def _cjk_font_files() -> tuple[str, ...]:
    """Font files from fc-list :lang=zh plus the usual XDG/system font dirs.

    Cached for the process lifetime: the scan rglobs the whole font tree and is
    otherwise repeated on every plot call.
    """
    import subprocess
    from pathlib import Path

    files: list[str] = []
    try:
        r = subprocess.run(["fc-list", ":lang=zh", "file"], capture_output=True, text=True, timeout=10)
        files += [l.strip().rstrip(":") for l in r.stdout.splitlines() if l.strip()]
    except Exception:
        pass
    for d in ("~/.local/share/fonts", "~/.fonts", "/usr/share/fonts", "/usr/local/share/fonts"):
        files += [str(x) for x in Path(d).expanduser().rglob("*") if x.suffix.lower() in (".ttf", ".otf", ".ttc")]
    return tuple(files)


def use_cjk_font() -> None:
    """Pick an installed font that renders both CJK and latin (no tofu boxes).

    Warns when no such font exists: matplotlib would then draw every CJK
    character as an empty box.
    """
    import sys

    import matplotlib
    from matplotlib import font_manager
    from matplotlib.ft2font import FT2Font

    matplotlib.use("Agg")
    chars = [ord(c) for c in _PROBE]

    def covers(path: str) -> bool:
        try:
            font = FT2Font(path)
            return all(font.get_char_index(c) for c in chars)
        except Exception:
            return False

    found = False
    for path in _cjk_font_files():
        if not covers(path):
            continue
        try:
            font_manager.fontManager.addfont(path)
        except Exception:
            continue
        fam = FT2Font(path).family_name
        if fam:
            import matplotlib.pyplot as plt

            plt.rcParams["font.family"] = "sans-serif"
            plt.rcParams["font.sans-serif"] = [fam, "DejaVu Sans"]
            plt.rcParams["axes.unicode_minus"] = False
            found = True
            return
    if not found:
        print(
            "[qenrich] warning: no installed font covers Chinese; CJK labels will "
            "render as boxes (install e.g. Noto Sans CJK / Source Han Sans)",
            file=sys.stderr,
        )
