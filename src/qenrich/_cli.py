"""qenrich CLI: gene enrichment (GO / KEGG / Pfam / InterPro) for non-model organisms."""

import argparse
import json
import re
import sys
from pathlib import Path
from zipfile import BadZipFile

import pandas as pd

from . import __version__
from ._genelist import _WEIGHTED_CELL
from ._io import CACHE_SCHEMA, _PROP_STEM, atomic_write, cache_dir_for, cache_fresh, load_net, \
    load_object, open_text, read_names, safe_names, save_net, save_objects
from ._parsers import PARSERS
from ._plot import plot_heatmap, plot_results
from ._plot_enrichplot import plot_results_enrichplot
from ._sniff import FORMAT_LABELS, FORMATS, sniff


class _Objects:
    """Feature lookup over a run's parsed/cached objects, loading on demand.

    A run needs one feature; a 30k-gene annotation caches five objects and reading
    them all costs ~0.7s. ``counts`` is what the run log reports, filled from the
    cache stamp when the objects are not in memory.
    """

    def __init__(self, counts: dict[str, dict]):
        self.counts = counts

    def load(self, feature: str) -> pd.DataFrame:
        raise NotImplementedError


class _FromCache(_Objects):
    def __init__(self, cdir: Path, names: list[str], counts: dict[str, dict]):
        # a stamp from before the counts key reports no objects; the log then
        # falls back to the loaded net's counts
        super().__init__(counts or {n: {} for n in names})
        self.cdir, self.names = cdir, names

    def load(self, feature: str) -> pd.DataFrame:
        if feature not in self.names:
            raise KeyError(f"feature '{feature}' not in parsed objects {sorted(self.names)}; "
                           f"use -f one of {sorted(self.names)}")
        return load_object(self.cdir, feature)


class _InMemory(_Objects):
    def __init__(self, objects: dict[str, pd.DataFrame]):
        super().__init__({k: {"terms": v["source"].nunique(), "genes": v["target"].nunique()}
                          for k, v in objects.items()})
        self.objects = objects

    def load(self, feature: str) -> pd.DataFrame:
        if feature not in self.objects:
            raise KeyError(f"feature '{feature}' not in parsed objects {sorted(self.objects)}; "
                           f"use -f one of {sorted(self.objects)}")
        return self.objects[feature]


def _resolve_input(args) -> tuple[str, _Objects, str, Path | None]:
    """Resolve -i into (default_feature, objects, kind, cache_dir).

    ``net`` means the input is a single net with no feature to choose: a net TSV
    file, or one object file resolved through ``--db``. The object name is kept
    as ``default_feature`` so it shows up in the run log and in the message that
    rejects a conflicting ``-f``. ``cache_dir`` is set when the input has one
    (annotation file), so the OBO and propagated-net caches can live beside it.
    """
    inp = args.i
    no_cache = args.no_cache
    if Path(inp).is_file():
        fmt = args.format or sniff(inp)
        print(f"[qenrich] input: {inp} (detected: {FORMAT_LABELS.get(fmt, fmt)})")
        if fmt == "net":
            return "net", _InMemory(PARSERS["net"](inp)), "net", None
        cdir = cache_dir_for(inp)
        objects = None
        if not no_cache and cache_fresh(cdir, inp, fmt) and not args.eggnog_lvl:
            names = _cache_names(cdir)
            if names:  # a stamp whose objects are all missing is no cache at all
                print(f"[qenrich] using cache: {cdir}")
                counts = _cache_counts(cdir, names)
                objects = _FromCache(cdir, names, counts)
        if objects is None:
            kwargs = {"annot_lvl": args.eggnog_lvl} if fmt == "eggnog" else {}
            parsed = PARSERS[fmt](inp, **kwargs)
            if not parsed:  # before save_objects: an empty parse must not poison the cache
                raise ValueError(f"no annotations found in {inp}")
            if args.eggnog_lvl:  # don't cache level-filtered results
                print("[qenrich] parsed (level-filtered, not cached)")
            elif no_cache:
                print("[qenrich] parsed (--no-cache: cache neither read nor written)")
            else:
                try:
                    save_objects(parsed, cdir, inp, fmt)
                    print(f"[qenrich] parsed and cached: {cdir}")
                except OSError as e:  # read-only input dir must not abort the run
                    print(f"[qenrich] warning: could not write cache {cdir} ({e}); continuing",
                          file=sys.stderr)
            objects = _InMemory(parsed)
        for k, c in objects.counts.items():
            if c.get("terms") is None:  # stamp without counts: name the object only
                print(f"[qenrich] object '{k}'")
            else:
                print(f"[qenrich] object '{k}': {c['terms']} terms, {c['genes']} genes")
        feature = "go" if "go" in objects.counts else next(iter(objects.counts))
        return feature, objects, "annot", cdir
    obj_file = Path(args.db) / f"{Path(inp).stem}.tsv"
    if obj_file.is_file():
        print(f"[qenrich] object file: {obj_file}")
        name = Path(inp).stem
        return name, _InMemory({name: PARSERS["net"](obj_file)["net"]}), "net", None
    raise FileNotFoundError(
        f"-i must be an existing file or an object name in --db (looked for {obj_file})")


_PROP_NPZ = "go_propagated.npz"
_PROP_META = "go_propagated.json"


def _cache_names(cdir: Path) -> list[str]:
    """Object names a fresh cache serves: the stamp's manifest, not the directory.

    TSVs a previous parse left behind (a --format switch, an annotation whose
    objects changed) must not stay selectable; names whose file went missing are
    dropped. A legacy stamp without a manifest falls back to the directory.
    """
    try:
        named = sorted(json.loads((cdir / "meta.json").read_text()).get("objects", {}))
    except (OSError, json.JSONDecodeError):
        named = []
    if not named:
        named = sorted(p.stem for p in cdir.glob("*.tsv") if p.stem != _PROP_STEM)
    return [n for n in named if (cdir / f"{n}.tsv").is_file()]


def _cache_counts(cdir: Path, names: list[str]) -> dict[str, dict]:
    """Per-object term/gene counts from the cache stamp, so the run log does not
    have to read every object.

    A cache written before the stamp held counts is read once and the stamp is
    updated in place; a read-only cache dir just keeps reporting net counts.
    """
    meta_path = cdir / "meta.json"
    try:
        counts = json.loads(meta_path.read_text()).get("objects", {})
    except (OSError, json.JSONDecodeError):
        return {}
    if counts:
        return counts
    counts = {}
    for n in names:
        df = load_object(cdir, n)
        counts[n] = {"terms": int(df["source"].nunique()), "genes": int(df["target"].nunique())}
    try:
        meta = json.loads(meta_path.read_text())
        meta["objects"] = counts
        meta_path.write_text(json.dumps(meta))
    except (OSError, json.JSONDecodeError):
        pass
    return counts


def _prop_cache_fresh(prop_path: Path, obo_path: str) -> bool:
    """True when the cached propagated net was built from this net, OBO and qenrich.

    The propagated net is a function of the parsed annotation, the OBO and the
    code, so the stamp carries the annotation cache's own stamp (source mtime,
    size, format and version -- a ``--format`` switch re-parses the same file)
    plus the cache schema, the OBO path and its mtime.
    """
    meta = prop_path.parent / _PROP_META
    if not meta.is_file():
        return False
    try:
        stamp = json.loads(meta.read_text())
        source = json.loads((prop_path.parent / "meta.json").read_text())
        return (stamp.get("schema") == CACHE_SCHEMA
                and stamp["obo"] == str(Path(obo_path).resolve())
                and stamp["obo_mtime"] == Path(obo_path).stat().st_mtime
                and stamp["version"] == __version__
                and stamp["source_mtime"] == source.get("mtime")
                and stamp["source_size"] == source.get("size")
                and stamp["source_format"] == source.get("format")
                and stamp["source_version"] == source.get("version"))
    except (KeyError, json.JSONDecodeError, OSError):
        return False


def _stamp_prop_cache(prop_path: Path, obo_path: str) -> None:
    source = json.loads((prop_path.parent / "meta.json").read_text())
    stamp = {
        "schema": CACHE_SCHEMA,
        "obo": str(Path(obo_path).resolve()),
        "obo_mtime": Path(obo_path).stat().st_mtime,
        "version": __version__,
        "source_mtime": source.get("mtime"),
        "source_size": source.get("size"),
        "source_format": source.get("format"),
        "source_version": source.get("version"),
    }
    atomic_write(prop_path.parent / _PROP_META, lambda p: p.write_text(json.dumps(stamp)))


def _bundled_obo() -> str | None:
    """The go-basic.obo shipped with the package (or a source checkout root)."""
    here = Path(__file__).resolve().parent
    root = here.parent.parent
    for cand in (here / "data" / "go-basic.obo.gz",
                 root / "data" / "go-basic.obo",
                 root / "data" / "go-basic.obo.gz"):
        if cand.is_file():
            return str(cand)
    return None


def _read_bg(path: str) -> list[str]:
    """Read a background gene list.

    ``gene,weight`` / ``gene;weight`` pairs keep only the id (same markup gene
    lists accept); bare comma/semicolon lists are split. The pair form is only
    taken when the id half holds no separator, so an all-numeric comma list such
    as ``7157,672,675,1234`` stays four genes instead of one.
    """
    vals = []
    with open_text(path) as fh:
        for line in fh:
            for tok in line.split():
                m = _WEIGHTED_CELL.match(tok)
                if m:
                    vals.append(m.group(1))
                else:
                    vals.extend(t for t in tok.replace(",", " ").replace(";", " ").split() if t)
    return vals


def _select_columns(header, sets, spec):
    """Filter gene-list columns to those named in ``spec`` (comma-separated
    header names or 1-based indices); default (None) keeps all."""
    if not spec:
        return sets
    picks = [p.strip() for p in spec.split(",")]
    selected = []
    for p in picks:
        if p in header:
            selected.append(p)
        elif p.isdigit() and 1 <= int(p) <= len(header):
            selected.append(header[int(p) - 1])
        else:
            raise ValueError(f"unknown column '{p}'; header is {header}")
    return {k: v for k, v in sets.items() if k in selected}


def _name_columns(df: pd.DataFrame, en_of, zh_of) -> pd.DataFrame:
    """Insert English ``name`` (always, when resolvable) and Chinese ``name_zh``
    (when the caller supplied Chinese names via --zh) after ``term``."""
    zhs = [zh_of(t) for t in df["term"]]
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
    from ._enrich import drop_parents, prune_es_wide, run_ora, strip_suffix
    from ._genelist import read_genelist
    from ._obo import GeneOntology

    default_feature, objects, kind, cdir = _resolve_input(args)
    if kind == "net":
        if args.feature and args.feature != default_feature:
            raise ValueError(
                f"-f {args.feature} conflicts with the input, which is the single '{default_feature}' "
                f"net; pass that object as -i, or give the full annotation file to pick a feature")
        feature = default_feature
    else:
        feature = args.feature or default_feature

    go = None
    if args.no_obo and args.obo:
        print(f"[qenrich] warning: --no-obo wins; --obo {args.obo} ignored", file=sys.stderr)
    if args.obo == "":
        print("[qenrich] warning: --obo is empty; using the bundled OBO", file=sys.stderr)
        args.obo = None
    obo_path = None if args.no_obo else (args.obo or _bundled_obo())
    cacheable = not args.no_cache and not args.eggnog_lvl
    prop_path = cdir / _PROP_NPZ if (cdir is not None and cacheable) else None
    obo_cdir = cdir if cacheable else None  # --no-cache writes nothing, OBO pickle included
    # propagation is deterministic in (net, OBO, qenrich), so a cached result can
    # stand in for the raw net entirely when the feature is the go object it was
    # built from
    prop_fresh = (prop_path is not None and prop_path.is_file() and obo_path is not None
                  and _prop_cache_fresh(prop_path, obo_path))
    net = None
    if prop_fresh and feature == "go":
        try:
            net = load_net(prop_path)
            go = GeneOntology.cached(obo_path, obo_cdir)
            print(f"[qenrich] using propagated GO net: {prop_path}")
        except (OSError, ValueError, KeyError, IndexError, EOFError, BadZipFile) as e:
            # corrupt or foreign npz: rebuild it instead of failing the run
            print(f"[qenrich] warning: {prop_path} unreadable ({e}); re-propagating",
                  file=sys.stderr)
            go, prop_fresh = None, False
    if net is None:
        net = objects.load(feature)
    counts = objects.counts.get(feature) or {}
    print(f"[qenrich] enriching feature: {feature} ("
          f"{counts.get('terms', net['source'].nunique())} terms, "
          f"{counts.get('genes', net['target'].nunique())} genes)")

    if obo_path and not prop_fresh:
        # applies whenever the active net holds GO ids, whatever the input route
        # (annotation file, parsed object db, or net TSV)
        if net["source"].astype(str).str.match(r"GO:\d{7}$").any():
            onto = GeneOntology.cached(obo_path, obo_cdir)
            propagated = onto.propagate(net)
            if propagated.empty:
                # every id missing from this OBO (annotation newer than the OBO):
                # propagating would yield nothing, so keep the raw annotations
                print(f"[qenrich] warning: none of the {net['source'].nunique()} terms are in "
                      f"{obo_path}; skipping propagation", file=sys.stderr)
            else:
                go = onto
                net = propagated
                if prop_path is not None:
                    try:
                        save_net(propagated, prop_path)
                        _stamp_prop_cache(prop_path, obo_path)
                    except OSError as e:
                        print(f"[qenrich] warning: could not write {prop_path} ({e}); continuing",
                              file=sys.stderr)
                src = "bundled" if not args.obo else "given"
                print(f"[qenrich] propagated GO DAG ({src} OBO: {Path(obo_path).name}): "
                      f"{net['source'].nunique()} terms, {net['target'].nunique()} genes")
        else:
            # only the user's explicit request is worth a warning; the bundled
            # default silently steps aside for non-GO features
            if args.obo:
                print("[qenrich] warning: --obo needs a go net, ignoring", file=sys.stderr)
    elif obo_path:
        go = GeneOntology.cached(obo_path, obo_cdir)  # term names for the report

    header, sets, weighted = read_genelist(args.genelist, no_header=args.no_header)
    if weighted:
        print(f"[qenrich] note: column(s) {', '.join(repr(c) for c in weighted)} hold gene,weight "
              f"pairs; the weights are ignored, the gene ids are used",
              file=sys.stderr)
    sets = _select_columns(header, sets, args.columns)
    if args.strip_suffix:
        sets, net = strip_suffix(sets, net)
        print("[qenrich] stripped .N version suffixes from gene ids")

    # name resolvers. English: the OBO is authoritative for GO ids; the --zh
    # table's col 2 fills the rest (KEGG, Pfam). Chinese only on request: --zh
    # pulls in the bundled go-zh table, and a --zh TABLE's col 3 extends or
    # overrides it.
    en_map: dict[str, str] = {}
    zh_map: dict[str, str] = {}
    names_df = args._names_df
    if names_df is not None and names_df.shape[1] >= 2:
        en_map = dict(zip(names_df.iloc[:, 0], names_df.iloc[:, 1].fillna("")))
    if names_df is not None and names_df.shape[1] >= 3:
        zh_map = dict(zip(names_df.iloc[:, 0], names_df.iloc[:, 2].fillna("")))
    zh_src = args._zh_path
    bundled = _bundled_go_zh() if args.zh else None
    if bundled and (not zh_src or Path(zh_src).resolve() != Path(bundled).resolve()):
        zdf = read_names(bundled)
        base = dict(zip(zdf.iloc[:, 0], zdf.iloc[:, 2].fillna("")))
        base.update(zh_map)  # the user's table wins over the bundled one
        zh_map = base

    def en_of(term: str) -> str:
        # the OBO is the authority for GO ids; the --zh table fills the gaps (KEGG, Pfam)
        return (go.name(term) if go else "") or en_map.get(term) or ""

    def zh_of(term: str) -> str:
        return zh_map.get(term, "")

    def label_of(term: str) -> str:
        return zh_of(term) or en_of(term) or term

    outdir = Path(args.out or "qenrich_results")
    outdir.mkdir(parents=True, exist_ok=True)
    summary = []

    def report(results, stats, tag):
        safe = safe_names(results)
        for name, df in results.items():
            s = stats[name]
            if s["n_input"] and s["n_hit"] < 0.5 * s["n_input"]:
                print(
                    f"[qenrich] WARNING set '{name}': only {s['n_hit']}/{s['n_input']} genes found in the "
                    f"annotation universe — check gene ID style (version suffixes? use --strip-suffix)",
                    file=sys.stderr,
                )
            d = _name_columns(df, en_of, zh_of)
            # only the ORA path counts terms that passed tmin but hold no query gene
            detail = f"{s['n_pruned']} below tmin" + (f", {s['n_empty']} with no query gene" if "n_empty" in s else "")
            if d.empty:
                print(f"[qenrich] set '{name}' ({tag}): no terms to report ({detail})")
                continue
            d.to_csv(outdir / f"{safe[name]}_{tag}.tsv", sep="\t", index=False)
            d = d.copy()
            d.insert(0, "set", name)
            summary.append(d)
            nsig = int((df["padj"] < args.padj).sum())
            print(
                f"[qenrich] set '{name}' ({tag}): {s['n_input']} genes in, {s['n_hit']} in universe, "
                f"{s['n_terms']} terms tested ({detail}), {nsig} with padj<{args.padj}"
            )

    results, es_wide, stats = {}, pd.DataFrame(), {}
    bg = None
    if args.bg:
        bg = _read_bg(args.bg)
        if args.strip_suffix:
            bg = [re.sub(r"\.\d+$", "", g) for g in bg]
    if sets:
        results, es_wide, stats = run_ora(
            net, sets, tmin=args.tmin, bg=bg, alternative=args.alternative
        )
        if args.drop_parents:
            if go:
                results = drop_parents(results, go, thr=args.padj)
                es_wide = prune_es_wide(es_wide, results)
                print(f"[qenrich] dropped parent terms with significant children (padj<{args.padj})")
            else:
                print("[qenrich] warning: --drop-parents needs GO terms from an OBO (not --no-obo), ignoring",
                      file=sys.stderr)
        report(results, stats, "enrichment")

    if summary:
        pd.concat(summary).to_csv(outdir / "summary.tsv", sep="\t", index=False)
        print(f"[qenrich] results written to {outdir}")
    if args.plot and summary:
        label_map = None
        if args.labels == "name" and (go or names_df is not None):
            label_map = {t: label_of(t) for t in net["source"].unique()}
        # --labels id keeps ids on the heatmap too (label_of falls back to id)
        heat_label = (lambda t: t) if args.labels == "id" else label_of
        if args.style == "enrichplot":
            plot_results_enrichplot(outdir, results, stats, label_map, tag="ep")
        else:
            plot_results(outdir, results, es_wide, label_map)
        plot_heatmap(pd.concat(summary), outdir, heat_label)
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


def _bundled_go_zh() -> str | None:
    """The Chinese name table: cwd copy, repo-root copy, then the packaged one."""
    here = Path(__file__).resolve().parent
    root = here.parent.parent
    for cand in (Path.cwd() / "go_zh.tsv",
                 root / "go_zh.tsv",
                 here / "data" / "go_zh.tsv.gz"):
        if cand.is_file():
            return str(cand)
    return None


def _resolve_table(path: str) -> str:
    """--zh given by bare name (e.g. go_zh.tsv): fall back to the packaged copy."""
    p = Path(path)
    if p.is_file():
        return str(p)
    data = Path(__file__).resolve().parent / "data"
    for cand in (data / f"{p.name}.gz", data / p.name):
        if cand.is_file():
            print(f"[qenrich] --zh {path}: using packaged {cand.name}")
            return str(cand)
    return str(p)  # missing everywhere: read_names reports it


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="qenrich", description=__doc__)
    ap.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")
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
    ep.add_argument("--no-cache", action="store_true",
                    help="do not read or write the <annotation>.qenrich/ cache: parse from scratch "
                         "and leave no cache behind")
    ep.add_argument("--alternative", choices=["greater", "less"], default="greater",
                    help="ORA alternative hypothesis: 'greater' (default) is P(X >= k), the one-sided "
                         "over-representation test clusterProfiler's enrichGO runs; 'less' is P(X <= k), "
                         "which tests depletion")
    ep.add_argument("--obo", help="go-basic.obo: propagate parents + term names "
                    "(default: the go-basic.obo bundled in data/)")
    ep.add_argument("--no-obo", action="store_true",
                    help="skip the OBO: no DAG propagation, term ids instead of names")
    ep.add_argument("--zh", nargs="?", const="go_zh.tsv", default=None, metavar="TABLE",
                    help="add Chinese names: a name_zh column and Chinese plot labels from the "
                         "bundled go-zh table; give a TSV (col 2 English, col 3 Chinese) to merge "
                         "a custom table for other id spaces (KEGG, Pfam) or override entries")
    ep.add_argument("--tmin", type=int, default=5, help="drop terms with fewer targets (default 5)")
    ep.add_argument("--bg", help="background gene list file (default: all annotated genes)")
    ep.add_argument("--padj", type=float, default=0.05,
                    help="padj cutoff for counting significant terms (and for --drop-parents); "
                         "outputs are not filtered (default 0.05)")
    ep.add_argument("--strip-suffix", action="store_true", help="drop .N version suffixes from gene ids (Gene01.1 -> Gene01)")
    ep.add_argument("--drop-parents", action="store_true", help="GO only: collapse parent terms with significant children")
    ep.add_argument("--eggnog-lvl", help="eggNOG only: filter rows by max_annot_lvl")
    ep.add_argument("-o", "--out", help="output dir (default ./qenrich_results)")
    ep.add_argument("--format", choices=FORMATS, help="force input format")
    ep.add_argument("--no-header", action="store_true", help="gene list has no header row")
    ep.add_argument("--labels", choices=["name", "id"], default="name",
                    help="plot label style: term name (Chinese when available, else English) "
                         "or bare term id")
    ep.add_argument("--plot", action="store_true", help="write barplot/dotplot/heatmap PNGs")
    ep.add_argument("--style", choices=["matplotlib", "enrichplot"], default="matplotlib",
                    help="plot style: plain matplotlib (default) or enrichplot (GeneRatio dotplot, Count barplot, heatplot)")
    ap.set_defaults(func=cmd_enrich)

    args = ap.parse_args(argv)
    if args.cmd is None and (not args.i or not args.genelist):
        ap.error("enrichment requires -i INPUT and --genelist FILE (or use the 'parse' subcommand)")
    try:
        args._zh_path = _resolve_table(args.zh) if args.zh else None
        args._names_df = read_names(args._zh_path) if args._zh_path else None
        return args.func(args)
    except (ValueError, KeyError, FileNotFoundError, AssertionError, OSError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
