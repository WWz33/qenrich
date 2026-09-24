"""Plotting helpers: dc.pl wrappers for per-set plots plus the summary heatmap."""

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd


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


def plot_results(outdir: str | Path, results: dict, es_wide: pd.DataFrame, label_map: dict | None) -> None:
    """Barplot + dotplot per set via decoupler's dc.pl (no hand-rolled plotting)."""
    import matplotlib

    matplotlib.use("Agg")
    import decoupler as dc

    use_cjk_font()  # term labels may be Chinese (--desc go_zh.tsv)

    outdir = Path(outdir)
    safe = {name: name.replace("/", "_") for name in results}
    labels = _unique_label_map(es_wide.columns, label_map) if label_map else None
    if labels:
        es_wide = es_wide.rename(columns=labels)
    lab_len = max((_disp_len(str(c)) for c in es_wide.columns), default=0)
    figsize = (max(6.5, 3.0 + 0.11 * lab_len), 3.5)
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
        if labels:
            long["label"] = long["term"].map(labels).fillna(long["term"])
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
