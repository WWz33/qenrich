"""Plotting helpers: dc.pl wrappers for per-set plots plus the summary heatmap."""

from pathlib import Path

import numpy as np
import pandas as pd


def plot_results(outdir: str | Path, results: dict, es_wide: pd.DataFrame, label_map: dict | None) -> None:
    """Barplot + dotplot per set via decoupler's dc.pl (no hand-rolled plotting)."""
    import matplotlib

    matplotlib.use("Agg")
    import decoupler as dc

    use_cjk_font()  # term labels may be Chinese (--desc go_zh.tsv)

    outdir = Path(outdir)
    safe = {name: name.replace("/", "_") for name in results}
    if label_map:
        es_wide = es_wide.rename(columns=label_map)
    lab_len = max((len(str(c)) for c in es_wide.columns), default=0)
    figsize = (max(5.5, 2.0 + 0.06 * lab_len), 3.5)
    for name, df in results.items():
        if df.empty:
            continue
        fname = safe[name]
        dc.pl.barplot(es_wide.fillna(0.0), name=name, top=15, save=str(outdir / f"{fname}_barplot.png"),
                      dpi=150, figsize=figsize)
        long = df.head(15).copy()
        score_col = "nes" if "nes" in df.columns else "log_or"
        # floor keeps non-significant terms visible, no cap needed with scale=0.15
        long["nlp"] = (-np.log10(long["padj"].clip(lower=2.22e-16))).clip(lower=0.3)
        if label_map:
            long["label"] = long["term"].map(label_map)
        else:
            long["label"] = long["term"]
        dc.pl.dotplot(long, x=score_col, y="label", c=score_col, s="nlp", top=15, scale=0.15,
                      save=str(outdir / f"{fname}_dotplot.png"), dpi=150, figsize=figsize)


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
    fig, ax = plt.subplots(figsize=(1.0 + 1.4 * mat.shape[1], 0.22 * len(mat) + 1.0), dpi=150)
    im = ax.imshow(mat.values, aspect="auto", cmap="Reds")
    ax.set_xticks(range(mat.shape[1]), mat.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(mat)), [name_of(t) or t for t in mat.index], fontsize=6)
    fig.colorbar(im, ax=ax, label="-log10(padj)", shrink=0.6)
    fig.savefig(outdir / "summary_heatmap.png", bbox_inches="tight")


_PROBE = "细胞过程Aa1"  # must cover CJK + latin + digits in ONE font


def _cjk_font_files() -> list[str]:
    """Font files from fc-list :lang=zh plus the usual XDG/system font dirs."""
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
    return files


def use_cjk_font() -> None:
    """Pick an installed font that renders both CJK and latin (no tofu boxes)."""
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
            return
