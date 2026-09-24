import numpy as np
import pandas as pd
import scipy.stats as sts
from pathlib import Path

from qenrich._enrich import run_ora
from qenrich._genelist import read_genelist
from qenrich._io import cache_dir_for, cache_fresh, load_objects, save_objects
from qenrich._obo import GeneOntology
from qenrich._parsers import PARSERS
from qenrich._sniff import sniff
from tests.fixtures import (
    B2G_ANNOT,
    B2G_TABULAR,
    EGGNOG,
    GAF,
    GFF3,
    GENERIC,
    IPRSCAN,
    KOFAM,
    NET,
    PANNZER,
    TRINOTATE,
    GENELIST,
    GENELIST_NOHEADER,
    OBO,
)


def wfile(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


# ----------------------------------------------------------------- sniff ----
def test_sniff_all_formats(tmp_path):
    cases = {
        "a.emapper.annotations.tsv": (EGGNOG, "eggnog"),
        "b.iprscan.tsv": (IPRSCAN, "iprscan"),
        "c.annot": (B2G_ANNOT, "b2g_annot"),
        "d.b2g.tsv": (B2G_TABULAR, "b2g_tabular"),
        "e.pannzer.out": (PANNZER, "pannzer"),
        "f.trinotate.tsv": (TRINOTATE, "trinotate"),
        "g.gaf": (GAF, "gaf"),
        "h.gff3": (GFF3, "gff3"),
        "i.kofam.txt": (KOFAM, "kofam"),
        "j.generic.txt": (GENERIC, "generic"),
        "k.net.tsv": (NET, "net"),
    }
    for name, (text, expected) in cases.items():
        assert sniff(wfile(tmp_path, name, text)) == expected, name


# --------------------------------------------------------------- parsers ----
def test_parse_eggnog(tmp_path):
    objs = PARSERS["eggnog"](wfile(tmp_path, "e.tsv", EGGNOG))
    assert set(objs) == {"go", "kegg", "pathway", "pfam", "cog"}
    go = objs["go"]
    assert set(go["target"]) == {"Gene01", "Gene02", "Gene03", "Gene05", "Gene06"}
    assert ("GO:0000001", "Gene01") in zip(go["source"], go["target"], strict=True)
    assert objs["pfam"]["source"].isin({"PF00001", "PF00002", "PF00003"}).all()  # version stripped
    assert objs["kegg"]["source"].isin({"K00001", "K00002", "K00003", "K00005"}).all()  # ko: prefix stripped


def test_parse_iprscan(tmp_path):
    objs = PARSERS["iprscan"](wfile(tmp_path, "i.tsv", IPRSCAN))
    assert set(objs) == {"go", "interpro", "pfam"}
    assert set(objs["go"]["target"]) == {"Gene01", "Gene02", "Gene03"}
    assert set(objs["interpro"]["source"]) == {"IPR000001", "IPR000002"}
    assert (objs["pfam"]["source"] == "PF00001").all()


def test_parse_gaf_drops_not(tmp_path):
    objs = PARSERS["gaf"](wfile(tmp_path, "g.gaf", GAF))
    assert ("GO:0000001", "Gene04") not in zip(objs["go"]["source"], objs["go"]["target"], strict=True)
    assert set(objs["go"]["target"]) == {"Gene01", "Gene02", "Gene03"}


def test_parse_gff3(tmp_path):
    objs = PARSERS["gff3"](wfile(tmp_path, "h.gff3", GFF3))
    # mRNA IPR000002 rolls up to its parent gene and "gene:" prefix is stripped
    assert set(objs["go"]["target"]) == {"Gene01", "Gene02", "Gene03"}
    assert set(objs["interpro"]["target"]) == {"Gene01"}
    assert set(objs["interpro"]["source"]) == {"IPR000001", "IPR000002"}


def test_parse_trinotate_merges_go_cols(tmp_path):
    objs = PARSERS["trinotate"](wfile(tmp_path, "t.tsv", TRINOTATE))
    go = objs["go"]
    assert set(zip(go["source"], go["target"], strict=True)) == {
        ("GO:0000001", "Gene01"),
        ("GO:0000002", "Gene01"),
        ("GO:0000001", "Gene02"),
        ("GO:0000003", "Gene02"),
    }


def test_parse_generic(tmp_path):
    objs = PARSERS["generic"](wfile(tmp_path, "g.txt", GENERIC))
    assert set(objs) == {"go", "kegg", "interpro"}
    assert set(objs["kegg"]["source"]) == {"K00001", "K00002"}


# ------------------------------------------------------------------- obo ----
def test_obo_propagation_and_meta(tmp_path):
    go = GeneOntology.from_obo(wfile(tmp_path, "go.obo", OBO))
    assert go.name("GO:0000002") == "child process"
    assert go.namespace("GO:0000003") == "cellular_component"
    assert "GO:0000001" in go.ancestors("GO:0000003")
    assert "GO:0000004" not in go._meta  # obsolete dropped
    assert go._alt["GO:0000099"] == "GO:0000001"


# -------------------------------------------------------------- genelist ----
def test_genelist_header_and_noheader(tmp_path):
    _, sets, num = read_genelist(wfile(tmp_path, "l.txt", GENELIST))
    assert set(sets) == {"up", "down"} and not num
    assert sets["up"] == ["Gene01", "Gene02", "Gene03"]
    _, sets, _ = read_genelist(wfile(tmp_path, "l2.txt", GENELIST_NOHEADER))
    assert set(sets) == {"Gene01", "Gene04"}  # all-unique ids: auto-detect must be forced with --no-header
    _, sets, _ = read_genelist(wfile(tmp_path, "l3.txt", GENELIST_NOHEADER), no_header=True)
    assert set(sets) == {"set1", "set2"}


# -------------------------------------------------------------------- io ----
def test_cache_roundtrip_and_fresh(tmp_path):
    import pandas as pd

    objs = {"go": pd.DataFrame({"source": ["GO:1"], "target": ["G1"]})}
    src = wfile(tmp_path, "src.tsv", "x")
    cdir = cache_dir_for(src)
    save_objects(objs, cdir, src, "eggnog")
    assert cache_fresh(cdir, src)
    assert list(load_objects(cdir)) == ["go"]
    Path(src).write_text("changed")
    assert not cache_fresh(cdir, src)


# ---------------------------------------------------------------- enrich ----
def _bh(pvals):
    import numpy as np

    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    adj = np.empty(n)
    run = 1.0
    for r in range(n - 1, -1, -1):
        run = min(run, p[order[r]] * n / (r + 1))
        adj[order[r]] = run
    return adj


def test_ora_matches_fisher_exact(tmp_path):
    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    _, sets, _ = read_genelist(wfile(tmp_path, "l.txt", GENELIST))
    results, es_wide, stats = run_ora(net, sets, tmin=3)
    up = results["up"]
    # universe=6 genes, three 3-gene terms (GO:1 {G1,G2,G3}, GO:2 {G1,G3,G4}, GO:3 {G1,G5,G6})
    assert stats["up"]["n_hit"] == 3
    row = up[up["term"] == "GO:0000003"].iloc[0]  # overlap with up={G1,G2,G3} is {G1}
    assert row["overlap"] == 1 and row["term_size"] == 3
    assert row["genes"] == "Gene01"
    # pvalue = two-sided Fisher exact on [[k, term-k], [n-k, neither]]
    assert abs(row["pvalue"] - sts.fisher_exact([[1, 2], [2, 1]])[1]) < 1e-12
    # padj = BH over the same two-sided Fisher p-values
    fisher = [
        sts.fisher_exact([[3, 0], [0, 3]])[1],  # GO:0000001, k=3
        sts.fisher_exact([[2, 1], [1, 2]])[1],  # GO:0000002, k=2
        sts.fisher_exact([[1, 2], [2, 1]])[1],  # GO:0000003, k=1
    ]
    for term, expect in zip(["GO:0000001", "GO:0000002", "GO:0000003"], _bh(fisher), strict=True):
        got = up[up["term"] == term].iloc[0]["padj"]
        assert abs(got - expect) < 1e-6, (term, got, expect)
    # log_or = Haldane-Anscombe corrected log odds ratio
    lor = np.log((1 + 0.5) * (1 + 0.5) / ((2 + 0.5) * (2 + 0.5)))
    assert abs(row["log_or"] - lor) < 1e-9
    assert up[up["term"] == "GO:0000001"].iloc[0]["log_or"] > 0  # enriched
    assert set(es_wide.columns) == {"GO:0000001", "GO:0000002", "GO:0000003"}


def test_ora_bg_and_tmin(tmp_path):
    import pytest

    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    bg = ["Gene01", "Gene02", "Gene03", "Gene04"]
    results, _, stats = run_ora(net, {"s": ["Gene01", "Gene02"]}, tmin=3, bg=bg)
    assert stats["s"]["n_hit"] == 2
    sizes = results["s"].set_index("term")["term_size"].to_dict()
    assert sizes["GO:0000001"] == 3 and sizes["GO:0000002"] == 3
    assert "GO:0000003" not in sizes  # only Gene01 left in bg -> below tmin=3
    # tmin prunes every term -> empty result table, not a crash
    empty, _, st = run_ora(net, {"s": ["Gene01"]}, tmin=10)
    assert empty["s"].empty and st["s"]["n_pruned"] == 3


# ------------------------------------------------------- new v0.2 features ----
GENELIST_NUMERIC = "up\tfc\nGene01\tGene01,2.0\nGene02\tGene03,-1.5\nGene03\tGene05,0.8\n"
GENELIST_CRLF = "up\r\nGene01\r\nGene02\r\nGene03\r\n"
DESC = "K00001\tpyruvate kinase\nK00002\thexokinase\n"


def test_crlf_genelist_and_gaf(tmp_path):
    from tests.fixtures import GAF as GAF_TXT

    _, sets, _ = read_genelist(wfile(tmp_path, "crlf.txt", GENELIST_CRLF))
    assert sets["up"] == ["Gene01", "Gene02", "Gene03"]
    objs = PARSERS["gaf"](wfile(tmp_path, "crlf.gaf", GAF_TXT.replace("\n", "\r\n")))
    assert set(objs["go"]["target"]) == {"Gene01", "Gene02", "Gene03"}


def test_gzip_transparent(tmp_path):
    import gzip as gz

    p = tmp_path / "n.tsv.gz"
    p.write_text(NET) if False else gz.open(p, "wt").write(NET)
    assert sniff(str(p)) == "net"
    objs = PARSERS["net"](str(p))
    assert len(objs["net"]) == 9


def test_strip_suffix(tmp_path):
    from qenrich._enrich import strip_suffix

    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    net["target"] = [f"{g}.1" for g in net["target"]]
    sets = {"s": ["Gene01.1", "Gene02"]}
    sets2, num2, net2 = strip_suffix(sets, {}, net)
    assert set(net2["target"]) == {"Gene01", "Gene02", "Gene03", "Gene04", "Gene05", "Gene06"}
    assert sets2["s"] == ["Gene01", "Gene02"]


def test_eggnog_annot_lvl(tmp_path):
    text = EGGNOG.replace("K\tE\n", "K\tE\n").replace(
        "#query\tseed_ortholog", "#query\tmax_annot_lvl\tseed_ortholog"
    )
    # inject lvl column values: build a tiny file directly instead
    lines = [
        "## x", "#query\tmax_annot_lvl\tGOs",
        "Gene01\t2\tGO:0000001", "Gene02\t2759\tGO:0000002", "Gene03\t2759\tGO:0000003",
    ]
    p = wfile(tmp_path, "lvl.emapper.tsv", "\n".join(lines) + "\n")
    objs = PARSERS["eggnog"](p, annot_lvl=2759)
    assert set(objs["go"]["target"]) == {"Gene02", "Gene03"}


def test_gsea_numeric_column(tmp_path):
    from qenrich._enrich import run_gsea

    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    _, sets, num = read_genelist(wfile(tmp_path, "num.txt", GENELIST_NUMERIC))
    assert set(sets) == {"up"} and set(num) == {"fc"}
    results, nes_wide, stats = run_gsea(net, num, tmin=1)  # 3-gene subset: any tmin>1 prunes all
    df = results["fc"]
    assert set(df.columns) == {"term", "term_size", "Count", "genes", "nes", "padj", "GeneRatio", "setSize"}
    assert len(df) == 3
    assert stats["fc"]["n_hit"] == 3
    assert df["Count"].notna().all() and (df["Count"] >= 0).all()
    assert df["nes"].notna().all()


def test_drop_parents(tmp_path):
    from qenrich._enrich import drop_parents
    from qenrich._obo import GeneOntology

    go = GeneOntology.from_obo(wfile(tmp_path, "go.obo", OBO))
    df = pd.DataFrame({
        "term": ["GO:0000001", "GO:0000002"],
        "padj": [0.01, 0.02],  # both significant; GO:1 is parent of GO:2
    })
    out = drop_parents({"s": df}, go, thr=0.05)["s"]
    assert list(out["term"]) == ["GO:0000002"]  # parent collapsed
    df2 = pd.DataFrame({"term": ["GO:0000001", "GO:0000002"], "padj": [0.01, 0.5]})
    out2 = drop_parents({"s": df2}, go, thr=0.05)["s"]
    assert list(out2["term"]) == ["GO:0000001", "GO:0000002"]  # child not significant: nothing to collapse


def test_read_names(tmp_path):
    from qenrich._io import read_names

    df = read_names(wfile(tmp_path, "desc.tsv", DESC))
    assert df.shape[1] == 2 and df.iloc[0, 0] == "K00001" and df.iloc[0, 1] == "pyruvate kinase"


def test_obo_children(tmp_path):
    from qenrich._obo import GeneOntology

    go = GeneOntology.from_obo(wfile(tmp_path, "go.obo", OBO))
    assert go.children("GO:0000001") == {"GO:0000002", "GO:0000003"}
    assert go.children("GO:0000003") == set()


def test_cli_end_to_end_numeric_and_flags(tmp_path, capsys):
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text(NET)
    (d / "gl.txt").write_text("up\tfc\nGene01\tGene01,2.0\nGene02\tGene03,-1.5\nGene03\tGene05,0.8\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"), "--tmin", "1",
               "-o", str(d / "out"), "--strip-suffix"])
    assert rc == 0
    assert (d / "out" / "up_enrichment.tsv").is_file()
    assert (d / "out" / "fc_gsea.tsv").is_file()
    assert (d / "out" / "summary.tsv").is_file()


EGGNOG_V2 = (
    "## emapper-2.1.13\n"
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\t"
    "Preferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tPFAMs\n"
    "Gene01\t-\t1e-5\t100\tAT6@ku\t2759\tK\tkinase\t-\tGO:0000001,GO:0000002\t2.7.1.1\tko:K00001\tko00010\tProfilin,PF00001.10\n"
    "Gene02\t-\t1e-9\t200\tAT6@ku\t2759\tE\tsynth\t-\tGO:0000001\t-\tko:K00002\tko00010\tArg_tRNA_synt_N\n"
)

EGGNOG_V3 = (
    "## emapper-3.0.0-beta6\n## applied filters:\n##   annot_evalue=0.001\n"
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\ttax_ceiling\tCOG_category\tPreferred_name\t"
    "GOs\tKEGG_ko\tKEGG_Pathway\tPFAMs\tannotation_confidence\n"
    "Gene01\t-\t1e-5\t100\tAT6@ku\t2759\tK\t-\tGO:0000003\tko:K00003\tko00020\tProfilin\t0.99\n"
    "Gene02\t-\t1e-9\t200\tAT6@ku\t2759\t-\t-\tGO:0000003\t-\t-\t-\t0.98\n"
)


def test_eggnog_v2_and_v3_headers(tmp_path):
    objs2 = PARSERS["eggnog"](wfile(tmp_path, "v2.tsv", EGGNOG_V2))
    assert set(objs2) == {"go", "kegg", "pathway", "pfam", "cog"}
    assert set(objs2["pfam"]["source"]) == {"Profilin", "PF00001", "Arg_tRNA_synt_N"}  # names + accessions
    objs3 = PARSERS["eggnog"](wfile(tmp_path, "v3.tsv", EGGNOG_V3))
    assert sniff(wfile(tmp_path, "v3s.tsv", EGGNOG_V3)) == "eggnog"
    assert set(objs3["go"]["target"]) == {"Gene01", "Gene02"}
    assert set(objs3["pfam"]["source"]) == {"Profilin"}
    # tax_ceiling filter (v3 renamed max_annot_lvl)
    objs3f = PARSERS["eggnog"](wfile(tmp_path, "v3f.tsv", EGGNOG_V3), annot_lvl=33090)
    assert not objs3f.get("go")
    objs3h = PARSERS["eggnog"](wfile(tmp_path, "v3h.tsv", EGGNOG_V3), annot_lvl=2759)
    assert set(objs3h["go"]["target"]) == {"Gene01", "Gene02"}


def test_enrichplot_style(tmp_path):
    from qenrich._plot_enrichplot import plot_results_enrichplot

    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    results, _, stats = run_ora(net, {"up": ["Gene01", "Gene02", "Gene03"]}, tmin=3)
    out = tmp_path / "ep"
    out.mkdir()
    plot_results_enrichplot(out, results, stats, None, top=3)
    for f in ("up_ep_dotplot.png", "up_ep_barplot.png", "ep_heatplot.png"):
        assert (out / f).is_file() and (out / f).stat().st_size > 1000
    d = results["up"]
    assert (d["overlap"] <= d["term_size"]).all()


def test_select_columns(tmp_path):
    from qenrich._cli import _select_columns
    header = ["up", "down", "fc"]
    sets = {"up": ["g1"], "down": ["g2"]}
    numeric = {"fc": {"g1": 1.0}}
    s2, n2 = _select_columns(header, sets, numeric, "up,fc")
    assert set(s2) == {"up"} and set(n2) == {"fc"}
    s3, n3 = _select_columns(header, sets, numeric, "2")
    assert set(s3) == {"down"} and not n3
    s4, n4 = _select_columns(header, sets, numeric, None)
    assert set(s4) == {"up", "down"} and set(n4) == {"fc"}


def test_name_columns_dual_and_no_zh(tmp_path):
    from qenrich._cli import _name_columns
    import pandas as pd
    df = pd.DataFrame({"term": ["GO:1", "GO:2"], "term_size": [3, 2]})
    en = {"GO:1": "stress response", "GO:2": "binding"}
    zh = {"GO:1": "应激响应", "GO:2": "结合"}
    d = _name_columns(df, lambda t: en.get(t, ""), lambda t: zh.get(t, ""))
    assert list(d.columns)[:3] == ["term", "name", "name_zh"]
    assert d.loc[0, "name_zh"] == "应激响应"
    d2 = _name_columns(df, lambda t: en.get(t, ""), lambda t: zh.get(t, ""), no_zh=True)
    assert "name_zh" not in d2.columns and "name" in d2.columns


# ---- review: direct parser tests for b2g_annot, b2g_tabular, pannzer, kofam ----
def test_parse_b2g_annot_direct(tmp_path):
    objs = PARSERS["b2g_annot"](wfile(tmp_path, "b.annot", B2G_ANNOT))
    assert set(zip(objs["go"]["source"], objs["go"]["target"])) == {
        ("GO:0000001", "Gene01"), ("GO:0000002", "Gene01"),
        ("GO:0000001", "Gene02"), ("GO:0000003", "Gene03")}
    assert set(zip(objs["interpro"]["source"], objs["interpro"]["target"])) == {
        ("IPR000001", "Gene01"), ("IPR000002", "Gene02"), ("IPR000001", "Gene03")}


def test_parse_b2g_tabular_multi_go(tmp_path):
    text = "Sequence Name\tGO IDs\nG1\tGO:0000001, GO:0000002\nG2\tGO:0000001\n"
    objs = PARSERS["b2g_tabular"](wfile(tmp_path, "b.tsv", text))
    go = objs["go"]
    assert set(zip(go["source"], go["target"])) == {
        ("GO:0000001", "G1"), ("GO:0000002", "G1"), ("GO:0000001", "G2")}


def test_parse_pannzer_multi_go(tmp_path):
    # inject a multi-value GO cell to test the finditer explosion
    text = "qpid\tgo_id\nGene01\tGO:0000001,GO:0000002\nGene02\tGO:0000001\n"
    objs = PARSERS["pannzer"](wfile(tmp_path, "p_mv.out", text))
    assert set(zip(objs["go"]["source"], objs["go"]["target"])) == {
        ("GO:0000001", "Gene01"), ("GO:0000002", "Gene01"), ("GO:0000001", "Gene02")}


def test_parse_kofam_direct(tmp_path):
    objs = PARSERS["kofam"](wfile(tmp_path, "k.out", KOFAM))
    kegg = objs["kegg"]
    assert ("K00001", "Gene01") in zip(kegg["source"], kegg["target"], strict=True)
    assert "-" not in set(kegg["source"])


# ---- review: GSEA leading edge numerical verification ----
def test_gsea_leading_edge_branches():
    """Asymmetric cases where positive and negative branches give DIFFERENT edges."""
    from qenrich._enrich import _leading_edge
    # positive peak dominates: hit g1 ranked first (weight 5), hit g6 ranked low (1)
    # positive edge = {g1}; negative edge would be {g6} — they differ
    vec = {"g1": 5.0, "g2": 4.0, "g3": -1.0, "g4": -2.0, "g5": -3.0, "g6": 1.0}
    count, edge = _leading_edge(vec, {"g1", "g6"})
    assert count == 1 and edge == ["g1"], (count, edge)
    # negative peak dominates: misses ranked high, hit g6 at the bottom
    # negative edge = {g6}; positive edge would be {g1, g6} — they differ
    vec = {"g1": 1.0, "g2": 5.0, "g3": 4.0, "g4": -1.0, "g5": -2.0, "g6": -3.0}
    count, edge = _leading_edge(vec, {"g1", "g6"})
    assert count == 1 and edge == ["g6"], (count, edge)


# ---- review: 3-column read_names (EN + ZH) ----
def test_read_names_3col(tmp_path):
    from qenrich._io import read_names
    text = "GO:0000001\tstress response\t应激响应\nGO:0000002\tbinding\t结合\n"
    df = read_names(wfile(tmp_path, "3col.tsv", text))
    assert df.shape[1] == 3
    assert df.iloc[0, 1] == "stress response"
    assert df.iloc[0, 2] == "应激响应"


# ---- review: empty gene set branch ----
def test_ora_empty_set(tmp_path):
    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    results, _, stats = run_ora(net, {"empty": []}, tmin=3)
    assert results["empty"].empty
    assert stats["empty"]["n_hit"] == 0


# ---- review: CLI e2e with --desc + --no-name-zh + -c ----
def test_cli_desc_no_name_zh_columns(tmp_path):
    from qenrich._cli import main
    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text(NET)
    (d / "gl.txt").write_text("up\tdown\nGene01\tGene04\nGene02\tGene05\nGene03\tGene06\n")
    (d / "zh.tsv").write_text("GO:0000001\tstress\t应激\nGO:0000002\tbinding\t结合\n"
                              "GO:0000003\tgrandchild\t孙节点\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"), "--tmin", "3",
               "-c", "up", "--desc", str(d / "zh.tsv"), "--no-name-zh",
               "-o", str(d / "out")])
    assert rc == 0
    import csv
    hdr = next(csv.reader(open(d / "out" / "up_enrichment.tsv"), delimiter="\t"))
    assert "name" in hdr and "name_zh" not in hdr  # no-name-zh drops Chinese


# ---- round 2: OBO propagate alt_id ----
def test_obo_propagate_alt_id(tmp_path):
    from qenrich._obo import GeneOntology
    import pandas as pd
    go = GeneOntology.from_obo(wfile(tmp_path, "go.obo", OBO))
    # GO:0000099 is an alt_id of GO:0000001; feeding it through propagate should
    # translate to primary and include GO:0000001's ancestors
    net = pd.DataFrame({"source": ["GO:0000099"], "target": ["Gene01"]})
    out = go.propagate(net)
    sources = set(out["source"])
    assert "GO:0000001" in sources  # alt translated to primary
    assert "GO:0000099" not in sources  # alt_id itself not in output
    # GO:0000001 is a root (no ancestors), so propagate yields only the primary
    assert sources == {"GO:0000001"}


# ---- round 2: GFF3 multi-Parent ----
def test_gff3_multi_parent(tmp_path):
    text = "##gff-version 3\n" \
           "chr1\tg\tgene\t1\t100\t.\t+\t.\tID=gene:Gene01;Ontology_term=GO:0000001\n" \
           "chr1\tg\tmRNA\t1\t100\t.\t+\t.\tID=mRNA1;Parent=gene:Gene01,gene:Gene99;Ontology_term=GO:0000002\n"
    objs = PARSERS["gff3"](wfile(tmp_path, "multi.gff3", text))
    # mRNA with comma-separated Parent should roll up to FIRST parent (Gene01), not "gene:Gene01,gene:Gene99"
    assert "Gene01" in set(objs["go"]["target"])
    assert "Gene99" not in set(objs["go"]["target"])  # second parent not used
    assert ("GO:0000002", "Gene01") in zip(objs["go"]["source"], objs["go"]["target"], strict=True)


# ---- round 2: read_names literal NA preservation ----
def test_read_names_literal_na(tmp_path):
    from qenrich._io import read_names
    text = "GO:0000001\tstress\t应激\nGO:0000002\tNA\t结合\n"
    df = read_names(wfile(tmp_path, "na.tsv", text))
    assert df.iloc[1, 1] == "NA"  # literal "NA" preserved, not dropped to NaN


# ---- round 4: numeric gene IDs survive cache reload ----
def test_cache_numeric_gene_ids(tmp_path):
    import pandas as pd
    objs = {"kegg": pd.DataFrame({"source": ["K1", "K2"], "target": ["00123", "456"]})}
    src = wfile(tmp_path, "anno.txt", "x")
    cdir = cache_dir_for(src)
    save_objects(objs, cdir, src, "kofam")
    loaded = load_objects(cdir)["kegg"]
    assert set(loaded["target"]) == {"00123", "456"}  # leading zero preserved


# ---- round 4: strip-suffix applies to --bg ----
def test_ora_bg_strip_suffix(tmp_path):
    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    net["target"] = [f"{g}.1" for g in net["target"]]
    bg = [f"{g}.1" for g in ["Gene01", "Gene02", "Gene03", "Gene04", "Gene05", "Gene06"]]
    sets = {"s": ["Gene01.1", "Gene02"]}
    from qenrich._enrich import strip_suffix, run_ora
    sets2, _, net2 = strip_suffix(sets, {}, net)
    results, _, stats = run_ora(net2, sets2, tmin=1, bg=[g[:-2] for g in bg])
    assert stats["s"]["n_hit"] == 2  # bg stripped: no overlap error, hits found


# ---- round 4: header with trailing tab (empty padded cell) ----
def test_genelist_trailing_tab_header(tmp_path):
    text = "SetA\tSetB\t\nGene01\tGene04\t\nGene02\tGene05\t\n"
    _, sets, _ = read_genelist(wfile(tmp_path, "tt.txt", text))
    assert set(sets) == {"SetA", "SetB"}  # header detected despite trailing tab
    assert sets["SetA"] == ["Gene01", "Gene02"]


# ---- round 4: OBO EOF stanza alt_id commit ----
def test_obo_eof_alt_id(tmp_path):
    from qenrich._obo import GeneOntology
    text = "[Term]\nid: GO:9a\nname: last term\nnamespace: molecular_function\nalt_id: GO:9b"  # no trailing blank line
    go = GeneOntology.from_obo(wfile(tmp_path, "eof.obo", text))
    assert go._alt.get("GO:9b") == "GO:9a"  # EOF commit registers alt_id


# ---- round 4: GFF3/B2G trailing tab tolerated ----
def test_gff3_trailing_tab(tmp_path):
    text = "##gff-version 3\nchr1\tg\tgene\t1\t100\t.\t+\t.\tID=gene:G1;Ontology_term=GO:0000001\t\n"
    objs = PARSERS["gff3"](wfile(tmp_path, "tt.gff3", text))
    assert ("GO:0000001", "G1") in zip(objs["go"]["source"], objs["go"]["target"], strict=True)


# ---- round 4: read_names literal values with BOM ----
def test_read_names_bom(tmp_path):
    from qenrich._io import read_names
    p = tmp_path / "bom.tsv"
    p.write_bytes("\ufeffGO:1\tstress\t应激\n".encode("utf-8"))
    df = read_names(str(p))
    assert df.iloc[0, 0] == "GO:1"  # BOM stripped, id intact


# ---- round 5: parse_generic keeps first row of header-less files ----
def test_generic_headerless_first_row_kept(tmp_path):
    text = "G1\tGO:0000001\nG2\tGO:0000002\n"
    objs = PARSERS["generic"](wfile(tmp_path, "nh.txt", text))
    assert set(objs["go"]["target"]) == {"G1", "G2"}  # G1 was previously dropped


# ---- round 5: single-row gene list is a header-only file ----
def test_genelist_single_row_header(tmp_path):
    _, sets, _ = read_genelist(wfile(tmp_path, "one.txt", "DE_up\tDE_down\n"))
    assert set(sets) == {"DE_up", "DE_down"} and sets["DE_up"] == []


# ---- round 5: strip_suffix dedups versioned net rows ----
def test_strip_suffix_net_dedup(tmp_path):
    from qenrich._enrich import strip_suffix
    net = pd.DataFrame({"source": ["GO:1", "GO:1"], "target": ["Gene01.1", "Gene01.2"]})
    _, _, net2 = strip_suffix({}, {}, net)
    assert not net2.duplicated(subset=["source", "target"]).any()
    assert set(net2["target"]) == {"Gene01"}


# ---- round 5: run_gsea survives fewer genes than tmin ----
def test_gsea_below_tmin_no_crash(tmp_path):
    from qenrich._enrich import run_gsea
    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    results, _, stats = run_gsea(net, {"tiny": {"Gene01": 2.0}}, tmin=5)
    assert results["tiny"].empty  # skipped, not crashed


# ---- round 5: iprscan versioned pfam accessions ----
def test_iprscan_versioned_pfam(tmp_path):
    text = ("G1\t0123456789abcdef0123456789abcdef\t150\tPfam\tPF00001.20\tKinase\t1\t100\t"
            "1e-5\tT\t20240101\tIPR000001\tDomain\tGO:0000001\t-\n")
    objs = PARSERS["iprscan"](wfile(tmp_path, "v.tsv", text))
    assert set(objs["pfam"]["source"]) == {"PF00001"}  # version stripped


# ---- round 5: gzipped OBO readable ----
def test_obo_gzipped(tmp_path):
    import gzip as gz
    from qenrich._obo import GeneOntology
    p = tmp_path / "t.obo.gz"
    with gz.open(p, "wt") as f:
        f.write(OBO)
    go = GeneOntology.from_obo(str(p))
    assert go.name("GO:0000001") == "parent process"


# ---- round 6: ragged rows must not shift pandas columns (index_col=False) ----
def test_eggnog_trailing_tab_no_column_shift(tmp_path):
    base = "#query\tGOs\nG1\tGO:0000001\t \nG2\tGO:0000002\t \n"
    objs = PARSERS["eggnog"](wfile(tmp_path, "tt.tsv", base))
    assert set(objs["go"]["target"]) == {"G1", "G2"}
    assert set(objs["go"]["source"]) == {"GO:0000001", "GO:0000002"}


def test_trinotate_trailing_tab_keeps_gene_ids(tmp_path):
    text = ("#gene_id\ttranscript_id\tgene_ontology_blast\n"
            "G1\tG1.t1\tGO:0004674^F\t \nG2\tG2.t1\tGO:0004177^F\t \n")
    objs = PARSERS["trinotate"](wfile(tmp_path, "tt.tsv", text))
    assert set(objs["go"]["target"]) == {"G1", "G2"}  # gene ids, not transcript ids


# ---- round 6: gff3 trailing tab+space must not drop the row ----
def test_gff3_trailing_tab_space(tmp_path):
    text = ("##gff-version 3\n"
            "chr1\tg\tgene\t1\t100\t.\t+\t.\tID=gene:G1;Ontology_term=GO:0000001\t \n")
    objs = PARSERS["gff3"](wfile(tmp_path, "ts.gff3", text))
    assert set(objs["go"]["target"]) == {"G1"}


# ---- round 6: pipe-separated identifiers everywhere ----
def test_gaf_pipe_separated_go(tmp_path):
    text = ("!gaf-version: 2.2\n"
            + "\t".join(["UniProtKB", "P12345", "P12345", "", "GO:0004674|GO:0005524",
                         "F:", "UniProtKB", "", "", "P", "20240101", "UniProtKB", "", "", ""]) + "\n")
    objs = PARSERS["gaf"](wfile(tmp_path, "p.gaf", text))
    assert set(objs["go"]["source"]) == {"GO:0004674", "GO:0005524"}


def test_generic_pipe_separated_ids(tmp_path):
    objs = PARSERS["generic"](wfile(tmp_path, "p.txt", "G1\tGO:0004177|GO:0004674\n"))
    assert set(objs["go"]["source"]) == {"GO:0004177", "GO:0004674"}


def test_trinotate_pipe_in_go_cell(tmp_path):
    text = "#gene_id\tgene_ontology_pfam\nG1\tGO:0000001|GO:0000002\n"
    objs = PARSERS["trinotate"](wfile(tmp_path, "p.tsv", text))
    assert set(objs["go"]["source"]) == {"GO:0000001", "GO:0000002"}


def test_kofam_pipe_joined_kos(tmp_path):
    objs = PARSERS["kofam"](wfile(tmp_path, "p.out", "G1\tK00001|K00002\nG2\tK00003\n"))
    assert set(objs["kegg"]["source"]) == {"K00001", "K00002", "K00003"}


def test_eggnog_uppercase_header(tmp_path):
    text = "#QUERY\tGOS\nG1\tGO:0000001\nG2\t-\n"
    objs = PARSERS["eggnog"](wfile(tmp_path, "u.tsv", text))
    assert set(objs["go"]["target"]) == {"G1"}


# ---- round 6: b2g .annot extra description column tolerated ----
def test_b2g_annot_fifth_column(tmp_path):
    text = "G1\tGO:0004674;GO:0005524\tInterPro\tIPR000719\tprotein kinase\nG2\tGO:0004177\tInterPro\tIPR001966\n"
    objs = PARSERS["b2g_annot"](wfile(tmp_path, "d.annot", text))
    assert set(objs["go"]["source"]) == {"GO:0004674", "GO:0005524", "GO:0004177"}
    assert set(objs["interpro"]["source"]) == {"IPR000719", "IPR001966"}
    assert "pfam" not in objs  # free-text 'PF00069' in col 5 must not fabricate pairs


# ---- round 6: ORA decoupler-free statistics ----
def test_ora_set_equals_universe(tmp_path):
    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    universe = sorted(set(net["target"]))
    results, _, _ = run_ora(net, {"all": universe}, tmin=3)
    assert not results["all"].empty  # previously crashed in decoupler's _runora


def test_ora_pvalue_padj_consistent(tmp_path):
    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    _, sets, _ = read_genelist(wfile(tmp_path, "l.txt", GENELIST))
    up = run_ora(net, sets, tmin=3)[0]["up"]
    sig = up[up["pvalue"] < 0.05]
    assert (sig["padj"] >= sig["pvalue"]).all()  # BH never shrinks p-values
    assert (up["padj"] <= 1.0).all() and (up["padj"] > 0).all()


def test_gsea_skip_path_no_nameerror(tmp_path):
    from qenrich._enrich import run_gsea
    # each term shares only 2 targets with the row -> below tmin=5 -> skip cleanly
    net = pd.DataFrame([(t, g) for t, gs in {"T1": [f"G{i}" for i in range(5)],
                                             "T2": [f"G{i}" for i in range(5, 10)]}
                        .items() for g in gs], columns=["source", "target"])
    results, _, stats = run_gsea(net, {"r": {g: 1.0 for g in ["G0", "G1", "G5", "G6"]}}, tmin=5)
    assert results["r"].empty and stats["r"]["n_input"] == 4


# ---- round 7: set names with / must not break plot filenames ----
def test_plot_slash_set_name(tmp_path):
    from qenrich._plot import plot_results

    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    results, es_wide, _ = run_ora(net, {"salt/control": ["Gene01", "Gene02", "Gene03"]}, tmin=3)
    out = tmp_path / "o"
    out.mkdir()
    plot_results(out, results, es_wide, None)
    assert (out / "salt_control_barplot.png").is_file()
    assert (out / "salt_control_dotplot.png").is_file()


def test_enrichplot_slash_set_name_and_tag(tmp_path):
    from qenrich._plot_enrichplot import plot_results_enrichplot

    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    results, _, stats = run_ora(net, {"a/b": ["Gene01", "Gene02", "Gene03"]}, tmin=3)
    out = tmp_path / "ep"
    out.mkdir()
    plot_results_enrichplot(out, results, stats, None, top=3)
    assert (out / "a_b_ep_dotplot.png").is_file()
    assert (out / "ep_heatplot.png").is_file()  # default tag


# ---- round 7: --obo works via the parse -> --db object route ----
def test_cli_obo_via_object_db(tmp_path, capsys):
    import subprocess, sys as _sys
    from qenrich._cli import main

    annot = wfile(tmp_path, "a.tsv", "#query\tGOs\nG1\tGO:0000001\nG2\tGO:0000001\nG3\tGO:0000001\n"
                                    "G4\tGO:0000003\nG5\tGO:0000003\nG6\tGO:0000003\n")
    gl = wfile(tmp_path, "g.txt", "s\nG1\nG2\nG3\nG4\n")
    db = tmp_path / "db"
    rc = main(["parse", str(annot), "-o", str(db)])
    assert rc == 0
    rc = main(["-i", "go", "--genelist", str(gl), "--db", str(db),
               "--obo", wfile(tmp_path, "t.obo", OBO), "--tmin", "1", "-o", str(tmp_path / "r")])
    assert rc == 0
    out = "\n".join(capsys.readouterr().out.splitlines())
    assert "propagated GO DAG" in out  # was silently skipped before
    res = (tmp_path / "r" / "s_enrichment.tsv").read_text()
    assert "name" in res.splitlines()[0]  # obo term names filled


# ---- round 7: ORA and GSEA heatplots coexist (no overwrite) ----
def test_enrichplot_ora_gsea_heatplot_coexist(tmp_path):
    from qenrich._plot_enrichplot import plot_results_enrichplot

    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    results, _, stats = run_ora(net, {"up": ["Gene01", "Gene02", "Gene03"]}, tmin=3)
    out = tmp_path / "ep"
    out.mkdir()
    plot_results_enrichplot(out, results, stats, None, top=3, tag="ep")
    plot_results_enrichplot(out, results, stats, None, top=3, tag="gsea_ep")
    assert (out / "ep_heatplot.png").is_file()
    assert (out / "gsea_ep_heatplot.png").is_file()


# ================= review round 8: verified-fix regression tests =================

def test_read_bg_weighted_pairs_and_comma_lists(tmp_path):
    """--bg keeps the id from gene,weight / gene;weight but still splits bare
    comma/semicolon lists, including all-numeric ones."""
    from qenrich._cli import _read_bg

    p = wfile(tmp_path, "bg.txt",
              "Gene01,3.2\nGene02;1.5\nGene03,Gene04\nGene05\n7157,672,675,1234\ng1;g2\n")
    assert _read_bg(p) == ["Gene01", "Gene02", "Gene03", "Gene04", "Gene05",
                           "7157", "672", "675", "1234", "g1", "g2"]


def test_num_cell_regex_rejects_separator_in_id():
    """The gene,weight pattern must not swallow a separator inside the id half."""
    from qenrich._genelist import _NUM_CELL

    assert _NUM_CELL.match("Gene01,3.2").group(1) == "Gene01"
    assert _NUM_CELL.match("Gene01;3.2").group(1) == "Gene01"
    assert _NUM_CELL.match("Gene01,-1.5e-3").group(1) == "Gene01"
    assert _NUM_CELL.match("7157,672,675,1234") is None
    assert _NUM_CELL.match("a,b,c") is None


def test_plot_unique_label_map_disambiguates():
    """Labels must be unique even when a display name collides with another term id."""
    from qenrich._plot import _unique_label_map

    m = _unique_label_map(["GO:1", "GO:2", "GO:3"], {"GO:1": "X", "GO:2": "X", "GO:3": "Y"})
    assert m["GO:1"] == "X" and m["GO:2"] == "GO:2" and m["GO:3"] == "Y"
    # a label that equals another term's id must still not collide
    m2 = _unique_label_map(["GO:1", "GO:2"], {"GO:1": "GO:2", "GO:2": "GO:2"})
    assert len(set(m2.values())) == 2 and set(m2) == {"GO:1", "GO:2"}
    m3 = _unique_label_map(["A", "B", "C"], {"A": "C", "B": "C", "C": "C"})
    assert len(set(m3.values())) == 3
    assert all(len({_unique_label_map(cols, lm)[c] for c in cols}) == len(cols)
               for cols, lm in [(["a", "b", "c"], {"a": "Z", "b": "Z", "c": "Z"}),
                                (["a", "b"], {"a": "b", "b": "b"})])


def test_obo_propagate_warns_on_dropped_terms(tmp_path, capsys):
    """Annotations to GO ids absent from the OBO are dropped with a warning."""
    from qenrich._obo import GeneOntology

    go = GeneOntology.from_obo(wfile(tmp_path, "go.obo", OBO))
    net = pd.DataFrame({"source": ["GO:0000002", "GO:9999999"], "target": ["g1", "g2"]})
    out = go.propagate(net)
    assert "GO:9999999" not in set(out["source"])
    assert "no entry in the OBO" in capsys.readouterr().err


def test_cjk_font_scan_cached():
    """_cjk_font_files must be memoized (no rglob on every plot call)."""
    from qenrich._plot import _cjk_font_files, use_cjk_font

    assert hasattr(_cjk_font_files, "cache_info")
    _cjk_font_files.cache_clear()
    use_cjk_font()
    use_cjk_font()
    assert _cjk_font_files.cache_info().misses == 1  # second call is a hit


def test_cli_version_uses_package_version(capsys):
    import pytest

    from qenrich import __version__
    from qenrich._cli import main

    with pytest.raises(SystemExit):
        main(["-V"])
    assert __version__ in capsys.readouterr().out


def test_cache_write_failure_does_not_abort(tmp_path, monkeypatch, capsys):
    """A read-only input dir (cache write fails) must warn, not abort the run."""
    from qenrich import _cli

    d = tmp_path / "run"
    d.mkdir()
    (d / "a.tsv").write_text("#query\tGOs\nG1\tGO:0000001\nG2\tGO:0000001\nG3\tGO:0000001\n")
    (d / "g.txt").write_text("s\nG1\nG2\nG3\n")

    def boom(*a, **k):
        raise PermissionError("read-only")

    monkeypatch.setattr(_cli, "save_objects", boom)
    rc = _cli.main(["-i", str(d / "a.tsv"), "--genelist", str(d / "g.txt"),
                    "--tmin", "1", "-o", str(d / "out")])
    assert rc == 0
    assert (d / "out" / "s_enrichment.tsv").is_file()
    assert "could not write cache" in capsys.readouterr().err


def test_gsea_padj_matches_decoupler_bh(tmp_path):
    """run_gsea must pass decoupler's BH-adjusted p-values through unchanged."""
    import decoupler as dc

    from qenrich._enrich import run_gsea

    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    num = {"fc": {"Gene01": 2.0, "Gene03": -1.5, "Gene05": 0.8}}
    results, _, _ = run_gsea(net, num, tmin=1)
    got = results["fc"].set_index("term")

    genes = sorted(num["fc"])
    row = pd.DataFrame([[num["fc"][g] for g in genes]], index=["fc"], columns=genes)
    es, pv = dc.mt.gsea(row, net, tmin=1, empty=False, verbose=False)

    for term in es.columns:
        exp_nes = float(es.loc["fc", term])
        exp_padj = max(float(pv.loc["fc", term]), np.finfo(float).eps)  # run_gsea clips at eps
        assert abs(got.loc[term, "nes"] - exp_nes) < 1e-9
        assert abs(got.loc[term, "padj"] - exp_padj) < 1e-9


def test_leading_edge_matches_decoupler_esrank():
    """_leading_edge must agree with decoupler's own ES peak on random vectors."""
    import numpy as np

    from decoupler.mt._gsea import _esrank

    from qenrich._enrich import _leading_edge

    rng = np.random.default_rng(0)
    for _ in range(200):
        n_all = int(rng.integers(6, 30))
        genes = [f"g{i}" for i in range(n_all)]
        vec = {g: float(rng.normal()) for g in genes}
        n_set = int(rng.integers(2, n_all))
        targets = set(rng.choice(genes, size=n_set, replace=False))

        count, edge = _leading_edge(vec, targets)

        order = sorted(vec, key=lambda g: -vec[g])
        row = np.ascontiguousarray([vec[g] for g in order], dtype=np.float64)
        set_msk = np.ascontiguousarray([g in targets for g in order])
        mx, j, _ = _esrank(row, np.arange(n_all), set_msk, 1.0 / (n_all - n_set))
        # decoupler's own peak index j defines the leading edge (clusterProfiler semantics)
        exp_edge = ({g for g in order[: j + 1] if g in targets} if mx >= 0
                    else {g for g in order[j + 1:] if g in targets})
        assert set(edge) == exp_edge, (vec, targets, j, mx)
        assert count == len(exp_edge)


# ================= round-2 review fixes =================

def test_read_names_header_row_skipped(tmp_path):
    """A --desc file with a header row must not inject bogus name mappings."""
    from qenrich._io import read_names

    text = "id\tname\tzh\nGO:0000001\tstress\t应激\nGO:0000002\tbinding\t结合\n"
    df = read_names(wfile(tmp_path, "hdr.tsv", text))
    assert set(df.iloc[:, 0]) == {"GO:0000001", "GO:0000002"}
    # header-less files are untouched
    df2 = read_names(wfile(tmp_path, "nohdr.tsv", "GO:0000001\tstress\t应激\n"))
    assert df2.iloc[0, 1] == "stress"


def test_leading_edge_all_zero_term_weights():
    """A term whose genes all carry weight 0 has no direction: empty edge."""
    from qenrich._enrich import _leading_edge

    vec = {"a": 0.0, "b": 0.0, "c": 3.0, "d": -2.0, "e": 1.0}
    assert _leading_edge(vec, {"a", "b"}) == (0, [])
    # mixed (one zero-weight hit among non-zero): the zero adds nothing but the
    # term still has direction, so an edge exists
    count, edge = _leading_edge({"a": 0.0, "b": 5.0, "c": -1.0, "d": -1.0}, {"a", "b"})
    assert count >= 1 and "b" in edge


def test_gsea_generatio_is_count_over_set_size(tmp_path):
    """GeneRatio follows clusterProfiler: Count (leading edge) / setSize."""
    from qenrich._enrich import run_gsea

    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    num = {"fc": {"Gene01": 2.0, "Gene02": 1.0, "Gene03": -1.0}}
    results, _, _ = run_gsea(net, num, tmin=1)
    d = results["fc"]
    assert "setSize" in d.columns
    assert (d["GeneRatio"] == d["Count"] / d["setSize"]).all()
    assert (d["setSize"] == 3).all()


def test_bg_ignored_for_gsea_warns(tmp_path, capsys):
    """--bg applies to ORA only; GSEA runs must warn about the ignore."""
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text(NET)
    (d / "gl.txt").write_text("fc\nGene01,2.0\nGene03,-1.5\nGene05,0.8\n")
    (d / "bg.txt").write_text("Gene01\nGene02\nGene03\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
               "--tmin", "1", "--bg", str(d / "bg.txt"), "-o", str(d / "out")])
    assert rc == 0
    assert "--bg applies to ORA only" in capsys.readouterr().err


def test_labels_id_applies_to_heatmap(tmp_path, capsys):
    """--labels id must keep term ids on the summary heatmap too."""
    from qenrich import _cli

    captured = {}
    orig = _cli.plot_heatmap

    def spy(summary, outdir, name_of):
        captured["label"] = name_of(summary["term"].iloc[0])
        return orig(summary, outdir, name_of)

    _cli.plot_heatmap = spy
    try:
        d = tmp_path / "run"
        d.mkdir()
        (d / "net.tsv").write_text(NET)
        (d / "gl.txt").write_text("up\nGene01\nGene02\nGene03\n")
        (d / "desc.tsv").write_text("GO:0000001\tstress\t应激\n")
        rc = _cli.main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
                        "--tmin", "1", "--desc", str(d / "desc.tsv"),
                        "--labels", "id", "-o", str(d / "out"), "--plot"])
        assert rc == 0
        assert captured["label"] == "GO:0000001"  # id, not "stress"
    finally:
        _cli.plot_heatmap = orig


# ================= plot label fixes: CJK width, layout, EN name default =================

def test_disp_len_cjk():
    """CJK chars count ~1.7x so figure width grows for Chinese labels."""
    from qenrich._plot import _disp_len
    from qenrich._plot_enrichplot import _disp_len as _dl2

    assert _disp_len("abcd") == 4.0
    assert _disp_len("对胁迫的响应") > 9.0
    assert _dl2("abcd") == 4.0  # both modules share the rule


def test_dotplot_axes_after_labels(tmp_path):
    """layout='tight' must place the axes to the RIGHT of the y labels."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import qenrich._plot_enrichplot as PE
    from qenrich._plot import use_cjk_font

    use_cjk_font()
    d = pd.DataFrame({
        "term": ["GO:1", "GO:2"], "Description": ["细胞周期检查点信号传导", "对胁迫的响应"],
        "padj": [1e-4, 1e-3], "Count": [3, 2], "GeneRatio": [0.3, 0.2],
    }).sort_values("GeneRatio", ascending=False).iloc[::-1]
    cap = {}
    orig = plt.Figure.savefig

    def spy(self, path, **kw):
        self.canvas.draw()
        ax = self.axes[0]
        r = self.canvas.get_renderer()
        need = max(ax.get_window_extent().x0 - t.get_window_extent(r).x0
                   for t in ax.get_yticklabels())
        cap["ok"] = ax.get_window_extent().x0 - need >= 0
        return orig(self, path, **kw)

    plt.Figure.savefig = spy
    try:
        PE._dotplot(d, tmp_path / "p.png", "s", (6.5, 3.8))
    finally:
        plt.Figure.savefig = orig
    assert cap["ok"]  # labels sit fully left of the axes


def test_en_default_uses_bundled_obo(tmp_path, capsys):
    """The bundled go-basic.obo is used by default: it propagates the DAG and
    supplies English term names, with no --obo flag."""
    from qenrich import _cli

    cap = {}
    orig_hm = _cli.plot_heatmap
    _cli.plot_heatmap = lambda s, o, n: (cap.update(lab=n(s["term"].iloc[0])),
                                         orig_hm(s, o, n))[1]
    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text("source\ttarget\nGO:0006950\tg1\nGO:0006950\tg2\nGO:0006950\tg3\n")
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    try:
        rc = _cli.main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
                        "--tmin", "1", "-o", str(d / "out"), "--plot"])
        assert rc == 0
        out = capsys.readouterr().out
        assert "propagated GO DAG (bundled OBO" in out
        assert cap["lab"] == "response to stress"  # GO:0006950, from the OBO
    finally:
        _cli.plot_heatmap = orig_hm


def test_no_obo_leaves_ids_as_labels(tmp_path, capsys):
    """--no-obo means no ontology: no propagation and no term names, just ids.
    go_zh.tsv is never picked up implicitly (--desc is the only way in)."""
    from qenrich import _cli

    cap = {}
    orig_hm = _cli.plot_heatmap
    _cli.plot_heatmap = lambda s, o, n: (cap.update(lab=n(s["term"].iloc[0])),
                                         orig_hm(s, o, n))[1]
    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text("source\ttarget\nGO:0006950\tg1\nGO:0006950\tg2\nGO:0006950\tg3\n")
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    try:
        rc = _cli.main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
                        "--tmin", "1", "--no-obo", "-o", str(d / "out"), "--plot"])
        assert rc == 0
        out = capsys.readouterr().out
        assert "propagated GO DAG" not in out
        assert "go_zh.tsv" not in out  # no implicit name source
        assert cap["lab"] == "GO:0006950"  # bare id
    finally:
        _cli.plot_heatmap = orig_hm


def test_desc_go_zh_resolves_to_packaged_copy(tmp_path, monkeypatch):
    """`--desc go_zh.tsv` works without a git checkout: the packaged copy is used."""
    from qenrich import _cli

    monkeypatch.chdir(tmp_path)  # no go_zh.tsv here
    assert _cli._resolve_desc("go_zh.tsv").endswith("go_zh.tsv.gz")
    df = _cli.read_names(_cli._resolve_desc("go_zh.tsv"))
    assert df.shape[1] == 3 and len(df) > 30000  # id, english, chinese


def test_desc_gives_chinese_and_english(tmp_path):
    """--desc supplies the Chinese column and names for ids the OBO lacks (KEGG)."""
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text("source\ttarget\nK00001\tg1\nK00001\tg2\nK00001\tg3\n")
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    (d / "ko.tsv").write_text("K00001\tpyruvate kinase\t丙酮酸激酶\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
               "--tmin", "1", "--desc", str(d / "ko.tsv"), "-o", str(d / "out")])
    assert rc == 0
    import csv
    with open(d / "out" / "s_enrichment.tsv") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert rows[0][:3] == ["term", "name", "name_zh"]
    assert rows[1][:3] == ["K00001", "pyruvate kinase", "丙酮酸激酶"]


# ================= bundled OBO: default on, obsolete-term redirection =================

def test_bundled_obo_is_found():
    from qenrich._cli import _bundled_obo
    p = _bundled_obo()
    assert p and Path(p).is_file()
    assert Path(p).name.startswith("go-basic")


def test_obsolete_term_maps_to_replacement(tmp_path):
    """A retired GO term with replaced_by keeps its annotations via the replacement."""
    from qenrich._obo import GeneOntology

    obo = (
        "[Term]\nid: GO:0000001\nname: live parent\nnamespace: biological_process\n\n"
        "[Term]\nid: GO:0000002\nname: retired\nnamespace: biological_process\n"
        "is_obsolete: true\nreplaced_by: GO:0000001\n\n"
        "[Term]\nid: GO:0000003\nname: retired no replacement\nnamespace: biological_process\n"
        "is_obsolete: true\n"
    )
    go = GeneOntology.from_obo(wfile(tmp_path, "r.obo", obo))
    assert go._alt["GO:0000002"] == "GO:0000001"  # redirected
    assert "GO:0000003" not in go._meta and "GO:0000003" not in go._alt  # unmappable, dropped
    net = pd.DataFrame({"source": ["GO:0000002", "GO:0000003"], "target": ["g1", "g2"]})
    out = go.propagate(net)
    assert set(zip(out["source"], out["target"])) == {("GO:0000001", "g1")}


def test_no_obo_flag_skips_propagation(tmp_path, capsys):
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text("source\ttarget\nGO:0006950\tg1\nGO:0006950\tg2\nGO:0006950\tg3\n")
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
               "--tmin", "1", "--no-obo", "-o", str(d / "out")])
    assert rc == 0
    assert "propagated GO DAG" not in capsys.readouterr().out


# ============ review-of-review fixes: revert guard, warnings, precedence ============

ALL_UNKNOWN_NET = ("source\ttarget\n"
                   "GO:9000001\tg1\nGO:9000001\tg2\nGO:9000002\tg3\n")


def test_empty_propagation_restores_raw_net(tmp_path, capsys):
    """When NO term maps into the OBO, enrichment must still run on the raw net."""
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text(ALL_UNKNOWN_NET)
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    (d / "tiny.obo").write_text(
        "[Term]\nid: GO:0000001\nname: live\nnamespace: biological_process\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
               "--tmin", "1", "--obo", str(d / "tiny.obo"), "-o", str(d / "out")])
    assert rc == 0
    assert "skipping propagation" in capsys.readouterr().err
    res = pd.read_csv(d / "out" / "s_enrichment.tsv", sep="\t")
    assert len(res) == 2 and set(res["term"]) == {"GO:9000001", "GO:9000002"}  # raw net used


def test_drop_warning_lists_user_ids_and_gene_loss(tmp_path, capsys):
    """The drop warning names the ids as they appear in the file, and reports
    genes that lose all GO annotations."""
    from qenrich._obo import GeneOntology

    obo = ("[Term]\nid: GO:0000001\nname: live\nnamespace: biological_process\n")
    go = GeneOntology.from_obo(wfile(tmp_path, "o.obo", obo))
    net = pd.DataFrame({"source": ["GO:0000001", "GO:0000001", "GO:9000003"],
                        "target": ["g1", "g2", "g4"]})  # g4 only via unknown term
    out = go.propagate(net)
    err = capsys.readouterr().err
    assert "GO:9000003" in err  # the user's id, not a translated one
    assert "1 gene(s) left with no GO annotation" in err
    assert set(out["target"]) == {"g1", "g2"}  # g4 dropped from the universe


def test_no_obo_and_obo_conflict_warns(tmp_path, capsys):
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text(NET)
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
               "--tmin", "1", "--no-obo", "--obo", str(d / "fake.obo"),
               "-o", str(d / "out")])
    assert rc == 0
    assert "--no-obo wins" in capsys.readouterr().err


def test_desc_missing_file_is_clean_error(tmp_path, capsys):
    """A --desc that exists nowhere prints `error:`, not a traceback."""
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text(NET)
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
               "--tmin", "1", "--desc", "nowhere.tsv", "-o", str(d / "out")])
    assert rc == 1
    assert "error:" in capsys.readouterr().err


def test_obo_name_beats_desc_for_go_ids(tmp_path):
    """English names for GO ids come from the OBO; --desc only fills the gaps."""
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text("source\ttarget\nGO:0006950\tg1\nGO:0006950\tg2\nGO:0006950\tg3\n")
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    # deliberately wrong English name for a GO id the OBO covers
    (d / "wrong.tsv").write_text("GO:0006950\tWRONG NAME\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
               "--tmin", "1", "--desc", str(d / "wrong.tsv"), "-o", str(d / "out")])
    assert rc == 0
    res = pd.read_csv(d / "out" / "s_enrichment.tsv", sep="\t")
    assert res.iloc[0]["name"] == "response to stress"  # OBO wins


def test_obsolete_alt_id_and_chain_redirect(tmp_path):
    """alt_ids of retired terms and replaced_by chains resolve to the live term."""
    from qenrich._obo import GeneOntology

    obo = (
        "[Term]\nid: GO:0000001\nname: live target\nnamespace: biological_process\n\n"
        "[Term]\nid: GO:0000002\nname: retired mid\nnamespace: biological_process\n"
        "is_obsolete: true\nreplaced_by: GO:0000001\nalt_id: GO:0000008\n\n"
        "[Term]\nid: GO:0000003\nname: retired head\nnamespace: biological_process\n"
        "is_obsolete: true\nreplaced_by: GO:0000002\n"
    )
    go = GeneOntology.from_obo(wfile(tmp_path, "c.obo", obo))
    assert go._alt["GO:0000003"] == "GO:0000001"  # chain followed
    assert go._alt["GO:0000008"] == "GO:0000001"  # alt_id of a retired term
    net = pd.DataFrame({"source": ["GO:0000003", "GO:0000008"], "target": ["g1", "g2"]})
    out = go.propagate(net)
    assert set(zip(out["source"], out["target"])) == {("GO:0000001", "g1"), ("GO:0000001", "g2")}
