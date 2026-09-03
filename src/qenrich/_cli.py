"""qenrich CLI: enrichment for non-model organisms, powered by decoupler."""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from ._io import cache_dir_for, cache_fresh, load_objects, open_text, read_names, save_objects
from ._parsers import PARSERS
from ._plot import plot_heatmap, plot_results
from ._plot_enrichplot import plot_results_enrichplot
from ._sniff import FORMAT_LABELS, FORMATS, sniff


def _resolve_input(args) -> tuple[str, dict[str, pd.DataFrame]]:
    """Resolve -i into (default_feature, objects); parse+cache raw annotation files."""
    inp = args.i
    if Path(inp).is_file():
        fmt = args.format or sniff(inp)
        print(f"[qenrich] input: {inp} (detected: {FORMAT_LABELS.get(fmt, fmt)})")
        if fmt == "net":
            return "net", PARSERS["net"](inp)
        cdir = cache_dir_for(inp)
        if cache_fresh(cdir, inp) and not getattr(args, "eggnog_lvl", None):
            print(f"[qenrich] using cache: {cdir}")
            objects = load_objects(cdir)
        else:
            kwargs = {"annot_lvl": args.eggnog_lvl} if fmt == "eggnog" else {}
            objects = PARSERS[fmt](inp, **kwargs)
            if not getattr(args, "eggnog_lvl", None):  # don't cache level-filtered results
                save_objects(objects, cdir, inp, fmt)
                print(f"[qenrich] parsed and cached: {cdir}")
            else:
                print(f"[qenrich] parsed (level-filtered, not cached)")
        if not objects:
            raise ValueError(f"no annotations found in {inp}")
        for k, v in objects.items():
            print(f"[qenrich] object '{k}': {v['source'].nunique()} terms, {v['target'].nunique()} genes")
        return ("go" if "go" in objects else next(iter(objects))), objects
    obj_file = Path(args.db) / f"{Path(inp).stem}.tsv"
    if obj_file.is_file():
        print(f"[qenrich] object file: {obj_file}")
        return "net", PARSERS["net"](obj_file)
    raise FileNotFoundError(
        f"-i must be an existing file or an object name in --db (looked for {obj_file})")


def _read_bg(path: str) -> list[str]:
    vals = []
    with open_text(path) as fh:
        for line in fh:
            vals.extend(t for t in line.strip().replace(",", " ").replace(";", " ").split() if t)
    return vals


def _select_columns(header, sets, numeric, spec):
    """Filter gene-list columns to those named in ``spec`` (comma-separated
    header names or 1-based indices); default (None) keeps all."""
    if not spec:
        return sets, numeric
    picks = [p.strip() for p in spec.split(",")]
    selected = []
    for p in picks:
        if p in header:
            selected.append(p)
        elif p.isdigit() and 1 <= int(p) <= len(header):
            selected.append(header[int(p) - 1])
        else:
            raise ValueError(f"unknown column '{p}'; header is {header}")
    return ({k: v for k, v in sets.items() if k in selected},
            {k: v for k, v in numeric.items() if k in selected})


def _name_columns(df: pd.DataFrame, en_of, zh_of, no_zh: bool = False) -> pd.DataFrame:
    """Insert English ``name`` (always, when resolvable) and Chinese ``name_zh``
    (unless ``no_zh`` / --no-name-zh) columns after ``term``."""
    zhs = [zh_of(t) for t in df["term"]] if not no_zh else []
    ens = [en_of(t) for t in df["term"]]
    if not any(zhs) and not any(ens):
        return df
    d = df.copy()
    if any(zhs):
        d.insert(1, "name_zh", zhs)
    if any(ens):
        d.insert(1, "name", ens)
    return d


def cmd_enrich(args) -> int:
    from ._enrich import drop_parents, run_gsea, run_ora, strip_suffix
    from ._genelist import read_genelist
    from ._obo import GeneOntology

    default_feature, objects = _resolve_input(args)
    if default_feature == "net":
        feature, net = "net", objects["net"]
    else:
        feature = args.feature or default_feature
        if feature not in objects:
            raise KeyError(f"feature '{feature}' not in parsed objects {sorted(objects)}; use -f one of {sorted(objects)}")
        net = objects[feature]
    print(f"[qenrich] enriching feature: {feature} ({net['source'].nunique()} terms, {net['target'].nunique()} genes)")

    go = None
    if args.obo:
        if feature != "go":
            print("[qenrich] warning: --obo only applies to the 'go' feature, ignoring", file=sys.stderr)
        else:
            go = GeneOntology.from_obo(args.obo)
            net = go.propagate(net)
            print(f"[qenrich] propagated GO DAG: {net['source'].nunique()} terms, {net['target'].nunique()} genes")

    header, sets, numeric = read_genelist(args.genelist, no_header=args.no_header)
    sets, numeric = _select_columns(header, sets, numeric, getattr(args, "columns", None))
    if args.strip_suffix:
        sets, numeric, net = strip_suffix(sets, numeric, net)
        print("[qenrich] stripped .N version suffixes from gene ids")

    # name resolvers: English from --desc col 2 or --obo; Chinese from --desc col 3
    names_df = args._names_df
    en_map, zh_map = {}, {}
    if names_df is not None:
        if names_df.shape[1] >= 2:
            en_map = dict(zip(names_df.iloc[:, 0], names_df.iloc[:, 1].fillna("")))
        if names_df.shape[1] >= 3:
            zh_map = dict(zip(names_df.iloc[:, 0], names_df.iloc[:, 2].fillna("")))

    def en_of(term: str) -> str:
        return en_map.get(term) or (go.name(term) if go else "") or ""

    def zh_of(term: str) -> str:
        return zh_map.get(term, "")

    def label_of(term: str) -> str:
        return zh_of(term) or en_of(term) or term

    outdir = Path(args.out or "qenrich_results")
    outdir.mkdir(parents=True, exist_ok=True)
    summary = []

    def report(results, stats, tag, score):
        for name, df in results.items():
            s = stats[name]
            if s["n_input"] and s["n_hit"] < 0.5 * s["n_input"]:
                print(
                    f"[qenrich] WARNING set '{name}': only {s['n_hit']}/{s['n_input']} genes found in the "
                    f"annotation universe — check gene ID style (version suffixes? use --strip-suffix)",
                    file=sys.stderr,
                )
            d = _name_columns(df, en_of, zh_of, no_zh=args.no_name_zh)
            if d.empty:
                print(f"[qenrich] set '{name}' ({tag}): no terms survived tmin={args.tmin}")
                continue
            safe = name.replace("/", "_")  # set names may contain / (header cells)
            d.to_csv(outdir / f"{safe}_{tag}.tsv", sep="\t", index=False)
            d = d.copy()
            d.insert(0, "set", name)
            summary.append(d)
            nsig = int((df["padj"] < args.padj).sum())
            print(
                f"[qenrich] set '{name}' ({tag}): {s['n_input']} genes in, {s['n_hit']} in universe, "
                f"{s['n_terms']} terms tested ({s['n_pruned']} pruned by tmin), {nsig} with padj<{args.padj}"
            )

    results, es_wide, stats = {}, pd.DataFrame(), {}
    results_gsea, nes_wide, stats_gsea = {}, pd.DataFrame(), {}
    bg = None
    if args.bg:
        bg = _read_bg(args.bg)
        if args.strip_suffix:
            import re as _re
            bg = [_re.sub(r"\.\d+$", "", g) for g in bg]
    if sets:
        results, es_wide, stats = run_ora(
            net, sets, tmin=args.tmin, bg=bg, verbose=args.verbose
        )
        if args.drop_parents:
            if go:
                results = drop_parents(results, go, thr=args.padj)
                print(f"[qenrich] dropped parent terms with significant children (padj<{args.padj})")
            else:
                print("[qenrich] warning: --drop-parents needs --obo, ignoring", file=sys.stderr)
        report(results, stats, "enrichment", "log_or")
    if numeric:
        results_gsea, nes_wide, stats_gsea = run_gsea(net, numeric, tmin=args.tmin, verbose=args.verbose)
        report(results_gsea, stats_gsea, "gsea", "nes")

    if summary:
        pd.concat(summary).to_csv(outdir / "summary.tsv", sep="\t", index=False)
        print(f"[qenrich] results written to {outdir}")
    if args.plot and summary:
        label_map = None
        if args.labels == "name" and (go or names_df is not None):
            label_map = {t: label_of(t) for t in net["source"].unique()}
        if args.style == "enrichplot":
            if sets:
                plot_results_enrichplot(outdir, results, stats, label_map)
            if numeric:
                plot_results_enrichplot(outdir, results_gsea, stats_gsea, label_map)
            plot_heatmap(pd.concat(summary), outdir, label_of)
        else:
            if sets:
                plot_results(outdir, results, es_wide, label_map)
            if numeric:
                plot_results(outdir, results_gsea, nes_wide, label_map)
            plot_heatmap(pd.concat(summary), outdir, label_of)
        print(f"[qenrich] plots written to {outdir}")
    return 0


def cmd_parse(args) -> int:
    fmt = args.format or sniff(args.annot)
    print(f"[qenrich] detected: {FORMAT_LABELS.get(fmt, fmt)}")
    kwargs = {"annot_lvl": args.eggnog_lvl} if fmt == "eggnog" else {}
    objects = PARSERS[fmt](args.annot, **kwargs)
    if not objects:
        raise ValueError(f"no annotations found in {args.annot}")
    outdir = Path(args.o or "qenrich_db")
    save_objects(objects, outdir, args.annot, fmt)
    for k, v in objects.items():
        print(f"[qenrich] object '{k}': {v['source'].nunique()} terms, {v['target'].nunique()} genes -> {outdir / (k + '.tsv')}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="qenrich", description=__doc__)
    ap.add_argument("-V", "--version", action="version", version="%(prog)s 0.1.0")
    sub = ap.add_subparsers(dest="cmd")

    pp = sub.add_parser("parse", help="parse an annotation file into object TSVs")
    pp.add_argument("annot", help="annotation file (auto-sniffed)")
    pp.add_argument("-o", help="output db dir (default ./qenrich_db)")
    pp.add_argument("--format", choices=FORMATS, help="force format")
    pp.add_argument("--eggnog-lvl", help="eggNOG only: filter rows by max_annot_lvl")
    pp.set_defaults(func=cmd_parse)

    ep = ap.add_argument_group("enrichment options")
    ep.add_argument("-i", help="annotation file, object name, or net TSV")
    ep.add_argument("--genelist", help="gene list file, one column per set")
    ep.add_argument("-c", "--columns", help="comma-separated gene-list columns to run (header names or 1-based indices; default: all)")
    ep.add_argument("-f", "--feature", help="object to enrich (go/kegg/pfam/interpro/cog/pathway); default: go if present")
    ep.add_argument("--db", default="qenrich_db", help="object db dir (default ./qenrich_db)")
    ep.add_argument("--obo", help="go-basic.obo: propagate parents + term names")
    ep.add_argument("--desc", help="id->names TSV; col 2 = English name, col 3 = Chinese name (e.g. go_zh.tsv)")
    ep.add_argument("--no-name-zh", action="store_true", help="omit the Chinese name column from output (English name stays)")
    ep.add_argument("--tmin", type=int, default=5, help="drop terms with fewer targets (default 5)")
    ep.add_argument("--bg", help="background gene list file (default: all annotated genes)")
    ep.add_argument("--padj", type=float, default=0.05, help="significance cutoff for summary (default 0.05)")
    ep.add_argument("--strip-suffix", action="store_true", help="drop .N version suffixes from gene ids (Gene01.1 -> Gene01)")
    ep.add_argument("--drop-parents", action="store_true", help="GO only: collapse parent terms with significant children")
    ep.add_argument("--eggnog-lvl", help="eggNOG only: filter rows by max_annot_lvl")
    ep.add_argument("-o", "--out", help="output dir (default ./qenrich_results)")
    ep.add_argument("--format", choices=FORMATS, help="force input format")
    ep.add_argument("--no-header", action="store_true", help="gene list has no header row")
    ep.add_argument("--labels", choices=["name", "id"], default="name",
                    help="plot label style: English term name (default; needs --obo/--desc, else falls back to id) or bare term id")
    ep.add_argument("--plot", action="store_true", help="write barplot/dotplot/heatmap PNGs")
    ep.add_argument("--style", choices=["matplotlib", "enrichplot"], default="matplotlib",
                    help="plot style: dc.pl matplotlib (default) or enrichplot (GeneRatio dotplot, Count barplot, heatplot)")
    ep.add_argument("-v", "--verbose", action="store_true")
    ap.set_defaults(func=cmd_enrich)

    args = ap.parse_args(argv)
    if args.cmd is None and (not args.i or not args.genelist):
        ap.error("enrichment requires -i INPUT and --genelist FILE (or use the 'parse' subcommand)")
    args._names_df = read_names(args.desc) if getattr(args, "desc", None) else None
    try:
        return args.func(args)
    except (ValueError, KeyError, FileNotFoundError, AssertionError, OSError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
