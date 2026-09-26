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
    _, sets, weighted = read_genelist(wfile(tmp_path, "l.txt", GENELIST))
    assert set(sets) == {"up", "down"} and not weighted
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
    # pvalue = phyper(k-1, M, N-M, n, lower.tail = FALSE), the one-sided tail
    assert abs(row["pvalue"] - sts.fisher_exact([[1, 2], [2, 1]], alternative="greater")[1]) < 1e-12
    # padj = BH over those p-values, for the terms holding at least one query gene
    fisher = [
        sts.fisher_exact([[3, 0], [0, 3]], alternative="greater")[1],  # GO:0000001, k=3
        sts.fisher_exact([[2, 1], [1, 2]], alternative="greater")[1],  # GO:0000002, k=2
        sts.fisher_exact([[1, 2], [2, 1]], alternative="greater")[1],  # GO:0000003, k=1
    ]
    for term, expect in zip(["GO:0000001", "GO:0000002", "GO:0000003"], _bh(fisher), strict=True):
        got = up[up["term"] == term].iloc[0]["padj"]
        assert abs(got - expect) < 1e-9, (term, got, expect)
    # log_or = Haldane-Anscombe corrected log odds ratio
    lor = np.log((1 + 0.5) * (1 + 0.5) / ((2 + 0.5) * (2 + 0.5)))
    assert abs(row["log_or"] - lor) < 1e-9
    assert up[up["term"] == "GO:0000001"].iloc[0]["log_or"] > 0  # enriched
    assert set(es_wide.columns) == {"GO:0000001", "GO:0000002", "GO:0000003"}


def test_ora_bg_and_tmin(tmp_path):
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


# --------------------------------------------------------------- features ----
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
    sets2, net2 = strip_suffix(sets, net)
    assert set(net2["target"]) == {"Gene01", "Gene02", "Gene03", "Gene04", "Gene05", "Gene06"}
    assert sets2["s"] == ["Gene01", "Gene02"]


def test_eggnog_annot_lvl(tmp_path):
    # build a tiny file with the level column directly
    lines = [
        "## x", "#query\tmax_annot_lvl\tGOs",
        "Gene01\t2\tGO:0000001", "Gene02\t2759\tGO:0000002", "Gene03\t2759\tGO:0000003",
    ]
    p = wfile(tmp_path, "lvl.emapper.tsv", "\n".join(lines) + "\n")
    objs = PARSERS["eggnog"](p, annot_lvl=2759)
    assert set(objs["go"]["target"]) == {"Gene02", "Gene03"}


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


def test_cli_end_to_end_weighted_column_and_flags(tmp_path, capsys):
    """A gene,weight column is accepted; the weights are dropped with a note."""
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text(NET)
    (d / "gl.txt").write_text("up\tweighted\nGene01\tGene01,2.0\nGene02\tGene03,-1.5\nGene03\tGene05,0.8\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"), "--tmin", "1",
               "-o", str(d / "out"), "--strip-suffix"])
    assert rc == 0
    assert (d / "out" / "up_enrichment.tsv").is_file()
    assert (d / "out" / "weighted_enrichment.tsv").is_file()
    assert (d / "out" / "summary.tsv").is_file()
    err = capsys.readouterr().err
    assert "weights are ignored" in err and "'weighted'" in err
    # the id half is what gets tested: every hit is one of those three ids, not
    # the "id,weight" text
    res = pd.read_csv(d / "out" / "weighted_enrichment.tsv", sep="\t")
    hit = {g for cell in res["genes"] for g in str(cell).split(";") if g}
    assert hit == {"Gene01", "Gene03", "Gene05"}


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
    header = ["up", "down", "third"]
    sets = {"up": ["g1"], "down": ["g2"], "third": ["g3"]}
    s2 = _select_columns(header, sets, "up,third")
    assert set(s2) == {"up", "third"}
    s3 = _select_columns(header, sets, "2")
    assert set(s3) == {"down"}
    s4 = _select_columns(header, sets, None)
    assert set(s4) == {"up", "down", "third"}


def test_name_columns(tmp_path):
    from qenrich._cli import _name_columns
    import pandas as pd
    df = pd.DataFrame({"term": ["GO:1", "GO:2"], "term_size": [3, 2]})
    en = {"GO:1": "stress response", "GO:2": "binding"}
    zh = {"GO:1": "应激响应", "GO:2": "结合"}
    d = _name_columns(df, lambda t: en.get(t, ""), lambda t: zh.get(t, ""))
    assert list(d.columns)[:3] == ["term", "name", "name_zh"]
    assert d.loc[0, "name_zh"] == "应激响应"
    # no Chinese names resolvable (no --zh): English only, no empty name_zh column
    d2 = _name_columns(df, lambda t: en.get(t, ""), lambda t: "")
    assert "name_zh" not in d2.columns and "name" in d2.columns


# ---- direct parser tests for b2g_annot, b2g_tabular, pannzer, kofam ----
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


# ---- 3-column read_names (EN + ZH) ----
def test_read_names_3col(tmp_path):
    from qenrich._io import read_names
    text = "GO:0000001\tstress response\t应激响应\nGO:0000002\tbinding\t结合\n"
    df = read_names(wfile(tmp_path, "3col.tsv", text))
    assert df.shape[1] == 3
    assert df.iloc[0, 1] == "stress response"
    assert df.iloc[0, 2] == "应激响应"


# ---- empty gene set branch ----
def test_ora_empty_set(tmp_path):
    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    results, _, stats = run_ora(net, {"empty": []}, tmin=3)
    assert results["empty"].empty
    assert stats["empty"]["n_hit"] == 0


# ---- CLI e2e with --zh + -c ----
def test_cli_desc_supplies_both_name_columns(tmp_path):
    from qenrich._cli import main
    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text(NET)
    (d / "gl.txt").write_text("up\tdown\nGene01\tGene04\nGene02\tGene05\nGene03\tGene06\n")
    (d / "zh.tsv").write_text("GO:0000001\tstress\t应激\nGO:0000002\tbinding\t结合\n"
                              "GO:0000003\tgrandchild\t孙节点\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"), "--tmin", "3",
               "-c", "up", "--zh", str(d / "zh.tsv"),
               "-o", str(d / "out")])
    assert rc == 0
    import csv
    hdr = next(csv.reader(open(d / "out" / "up_enrichment.tsv"), delimiter="\t"))
    assert "name" in hdr and "name_zh" in hdr  # --zh table col 3 brings Chinese


# ---- OBO propagate alt_id ----
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


# ---- GFF3 multi-Parent ----
def test_gff3_multi_parent(tmp_path):
    text = "##gff-version 3\n" \
           "chr1\tg\tgene\t1\t100\t.\t+\t.\tID=gene:Gene01;Ontology_term=GO:0000001\n" \
           "chr1\tg\tmRNA\t1\t100\t.\t+\t.\tID=mRNA1;Parent=gene:Gene01,gene:Gene99;Ontology_term=GO:0000002\n"
    objs = PARSERS["gff3"](wfile(tmp_path, "multi.gff3", text))
    # mRNA with comma-separated Parent should roll up to FIRST parent (Gene01), not "gene:Gene01,gene:Gene99"
    assert "Gene01" in set(objs["go"]["target"])
    assert "Gene99" not in set(objs["go"]["target"])  # second parent not used
    assert ("GO:0000002", "Gene01") in zip(objs["go"]["source"], objs["go"]["target"], strict=True)


# ---- read_names literal NA preservation ----
def test_read_names_literal_na(tmp_path):
    from qenrich._io import read_names
    text = "GO:0000001\tstress\t应激\nGO:0000002\tNA\t结合\n"
    df = read_names(wfile(tmp_path, "na.tsv", text))
    assert df.iloc[1, 1] == "NA"  # literal "NA" preserved, not dropped to NaN


# ---- numeric gene IDs survive cache reload ----
def test_cache_numeric_gene_ids(tmp_path):
    import pandas as pd
    objs = {"kegg": pd.DataFrame({"source": ["K1", "K2"], "target": ["00123", "456"]})}
    src = wfile(tmp_path, "anno.txt", "x")
    cdir = cache_dir_for(src)
    save_objects(objs, cdir, src, "kofam")
    loaded = load_objects(cdir)["kegg"]
    assert set(loaded["target"]) == {"00123", "456"}  # leading zero preserved


# ---- strip-suffix applies to --bg ----
def test_ora_bg_strip_suffix(tmp_path):
    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    net["target"] = [f"{g}.1" for g in net["target"]]
    bg = [f"{g}.1" for g in ["Gene01", "Gene02", "Gene03", "Gene04", "Gene05", "Gene06"]]
    sets = {"s": ["Gene01.1", "Gene02"]}
    from qenrich._enrich import strip_suffix, run_ora
    sets2, net2 = strip_suffix(sets, net)
    results, _, stats = run_ora(net2, sets2, tmin=1, bg=[g[:-2] for g in bg])
    assert stats["s"]["n_hit"] == 2  # bg stripped: no overlap error, hits found


# ---- header with trailing tab (empty padded cell) ----
def test_genelist_trailing_tab_header(tmp_path):
    text = "SetA\tSetB\t\nGene01\tGene04\t\nGene02\tGene05\t\n"
    _, sets, _ = read_genelist(wfile(tmp_path, "tt.txt", text))
    assert set(sets) == {"SetA", "SetB"}  # header detected despite trailing tab
    assert sets["SetA"] == ["Gene01", "Gene02"]


# ---- OBO EOF stanza alt_id commit ----
def test_obo_eof_alt_id(tmp_path):
    from qenrich._obo import GeneOntology
    text = "[Term]\nid: GO:9a\nname: last term\nnamespace: molecular_function\nalt_id: GO:9b"  # no trailing blank line
    go = GeneOntology.from_obo(wfile(tmp_path, "eof.obo", text))
    assert go._alt.get("GO:9b") == "GO:9a"  # EOF commit registers alt_id


# ---- GFF3/B2G trailing tab tolerated ----
def test_gff3_trailing_tab(tmp_path):
    text = "##gff-version 3\nchr1\tg\tgene\t1\t100\t.\t+\t.\tID=gene:G1;Ontology_term=GO:0000001\t\n"
    objs = PARSERS["gff3"](wfile(tmp_path, "tt.gff3", text))
    assert ("GO:0000001", "G1") in zip(objs["go"]["source"], objs["go"]["target"], strict=True)


# ---- read_names literal values with BOM ----
def test_read_names_bom(tmp_path):
    from qenrich._io import read_names
    p = tmp_path / "bom.tsv"
    p.write_bytes("\ufeffGO:1\tstress\t应激\n".encode("utf-8"))
    df = read_names(str(p))
    assert df.iloc[0, 0] == "GO:1"  # BOM stripped, id intact


# ---- parse_generic keeps first row of header-less files ----
def test_generic_headerless_first_row_kept(tmp_path):
    text = "G1\tGO:0000001\nG2\tGO:0000002\n"
    objs = PARSERS["generic"](wfile(tmp_path, "nh.txt", text))
    assert set(objs["go"]["target"]) == {"G1", "G2"}  # the header heuristic must keep a real first-row gene


# ---- single-row gene list is a header-only file ----
def test_genelist_single_row_header(tmp_path):
    _, sets, _ = read_genelist(wfile(tmp_path, "one.txt", "DE_up\tDE_down\n"))
    assert set(sets) == {"DE_up", "DE_down"} and sets["DE_up"] == []


# ---- strip_suffix dedups versioned net rows ----
def test_strip_suffix_net_dedup(tmp_path):
    from qenrich._enrich import strip_suffix
    net = pd.DataFrame({"source": ["GO:1", "GO:1"], "target": ["Gene01.1", "Gene01.2"]})
    _, net2 = strip_suffix({}, net)
    assert not net2.duplicated(subset=["source", "target"]).any()
    assert set(net2["target"]) == {"Gene01"}


# ---- iprscan versioned pfam accessions ----
def test_iprscan_versioned_pfam(tmp_path):
    text = ("G1\t0123456789abcdef0123456789abcdef\t150\tPfam\tPF00001.20\tKinase\t1\t100\t"
            "1e-5\tT\t20240101\tIPR000001\tDomain\tGO:0000001\t-\n")
    objs = PARSERS["iprscan"](wfile(tmp_path, "v.tsv", text))
    assert set(objs["pfam"]["source"]) == {"PF00001"}  # version stripped


# ---- gzipped OBO readable ----
def test_obo_gzipped(tmp_path):
    import gzip as gz
    from qenrich._obo import GeneOntology
    p = tmp_path / "t.obo.gz"
    with gz.open(p, "wt") as f:
        f.write(OBO)
    go = GeneOntology.from_obo(str(p))
    assert go.name("GO:0000001") == "parent process"


# ---- ragged rows must not shift pandas columns (index_col=False) ----
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


# ---- gff3 trailing tab+space must not drop the row ----
def test_gff3_trailing_tab_space(tmp_path):
    text = ("##gff-version 3\n"
            "chr1\tg\tgene\t1\t100\t.\t+\t.\tID=gene:G1;Ontology_term=GO:0000001\t \n")
    objs = PARSERS["gff3"](wfile(tmp_path, "ts.gff3", text))
    assert set(objs["go"]["target"]) == {"G1"}


# ---- pipe-separated identifiers everywhere ----
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


# ---- b2g .annot extra description column tolerated ----
def test_b2g_annot_fifth_column(tmp_path):
    text = "G1\tGO:0004674;GO:0005524\tInterPro\tIPR000719\tprotein kinase\nG2\tGO:0004177\tInterPro\tIPR001966\n"
    objs = PARSERS["b2g_annot"](wfile(tmp_path, "d.annot", text))
    assert set(objs["go"]["source"]) == {"GO:0004674", "GO:0005524", "GO:0004177"}
    assert set(objs["interpro"]["source"]) == {"IPR000719", "IPR001966"}
    assert "pfam" not in objs  # free-text 'PF00069' in col 5 must not fabricate pairs


# ---- ORA statistics ----
def test_ora_set_equals_universe(tmp_path):
    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    universe = sorted(set(net["target"]))
    results, _, _ = run_ora(net, {"all": universe}, tmin=3)
    assert not results["all"].empty  # a universe-sized set still gets a table


def test_ora_pvalue_padj_consistent(tmp_path):
    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    _, sets, _ = read_genelist(wfile(tmp_path, "l.txt", GENELIST))
    up = run_ora(net, sets, tmin=3)[0]["up"]
    sig = up[up["pvalue"] < 0.05]
    assert (sig["padj"] >= sig["pvalue"]).all()  # BH never shrinks p-values
    assert (up["padj"] <= 1.0).all() and (up["padj"] > 0).all()


# ---- set names with / must not break plot filenames ----
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


# ---- --obo works via the parse -> --db object route ----
def test_cli_obo_via_object_db(tmp_path, capsys):
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
    assert "propagated GO DAG" in out  # the object-db route must propagate too
    res = (tmp_path / "r" / "s_enrichment.tsv").read_text()
    assert "name" in res.splitlines()[0]  # obo term names filled


# ===================== regression tests =====================

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
    from qenrich._genelist import _WEIGHTED_CELL

    assert _WEIGHTED_CELL.match("Gene01,3.2").group(1) == "Gene01"
    assert _WEIGHTED_CELL.match("Gene01;3.2").group(1) == "Gene01"
    assert _WEIGHTED_CELL.match("Gene01,-1.5e-3").group(1) == "Gene01"
    assert _WEIGHTED_CELL.match("7157,672,675,1234") is None
    assert _WEIGHTED_CELL.match("a,b,c") is None


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


# ========== read_names header skip, --labels id on the heatmap ==========

def test_read_names_header_row_skipped(tmp_path):
    """A --zh table with a header row must not inject bogus name mappings."""
    from qenrich._io import read_names

    text = "id\tname\tzh\nGO:0000001\tstress\t应激\nGO:0000002\tbinding\t结合\n"
    df = read_names(wfile(tmp_path, "hdr.tsv", text))
    assert set(df.iloc[:, 0]) == {"GO:0000001", "GO:0000002"}
    # header-less files are untouched
    df2 = read_names(wfile(tmp_path, "nohdr.tsv", "GO:0000001\tstress\t应激\n"))
    assert df2.iloc[0, 1] == "stress"


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
                        "--tmin", "1", "--zh", str(d / "desc.tsv"),
                        "--labels", "id", "-o", str(d / "out"), "--plot"])
        assert rc == 0
        assert captured["label"] == "GO:0000001"  # id, not "stress"
    finally:
        _cli.plot_heatmap = orig


# ================= plot label fixes: CJK width, layout, EN name default =================

def test_disp_len_cjk():
    """CJK chars count ~1.7x so figure width grows for Chinese labels."""
    from qenrich._plot import _disp_len

    assert _disp_len("abcd") == 4.0
    assert _disp_len("对胁迫的响应") > 9.0


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
        # label_of prefers Chinese: "response to stress" (OBO) + 对胁迫的响应 (built-in)
        assert cap["lab"] == "response to stress"  # default: English from the OBO
    finally:
        _cli.plot_heatmap = orig_hm


def test_zh_flag_adds_chinese_labels(tmp_path, capsys):
    """Chinese appears only on request: --zh pulls in the bundled table."""
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
                        "--tmin", "1", "--zh", "-o", str(d / "out"), "--plot"])
        assert rc == 0
        assert cap["lab"] == "对胁迫的响应"  # bundled go_zh.tsv
        import csv
        with open(d / "out" / "s_enrichment.tsv") as fh:
            hdr = next(csv.reader(fh, delimiter="\t"))
        assert hdr[:3] == ["term", "name", "name_zh"]
    finally:
        _cli.plot_heatmap = orig_hm


def test_zh_table_resolves_to_packaged_copy(tmp_path, monkeypatch):
    """`--zh go_zh.tsv` works without a git checkout: the packaged copy is used."""
    from qenrich import _cli

    monkeypatch.chdir(tmp_path)  # no go_zh.tsv here
    assert _cli._resolve_table("go_zh.tsv").endswith("go_zh.tsv.gz")
    df = _cli.read_names(_cli._resolve_table("go_zh.tsv"))
    assert df.shape[1] == 3 and len(df) > 30000  # id, english, chinese


def test_zh_table_gives_chinese_and_english(tmp_path):
    """--zh TABLE supplies the Chinese column and names for ids the OBO lacks (KEGG)."""
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text("source\ttarget\nK00001\tg1\nK00001\tg2\nK00001\tg3\n")
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    (d / "ko.tsv").write_text("K00001\tpyruvate kinase\t丙酮酸激酶\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
               "--tmin", "1", "--zh", str(d / "ko.tsv"), "-o", str(d / "out")])
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


# ============ propagation revert guard, warnings, flag precedence ============

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


def test_zh_missing_file_is_clean_error(tmp_path, capsys):
    """A --zh table that exists nowhere prints `error:`, not a traceback."""
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text(NET)
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
               "--tmin", "1", "--zh", "nowhere.tsv", "-o", str(d / "out")])
    assert rc == 1
    assert "error:" in capsys.readouterr().err


def test_obo_name_beats_desc_for_go_ids(tmp_path):
    """English names for GO ids come from the OBO; the --zh table only fills the gaps."""
    from qenrich._cli import main

    d = tmp_path / "run"
    d.mkdir()
    (d / "net.tsv").write_text("source\ttarget\nGO:0006950\tg1\nGO:0006950\tg2\nGO:0006950\tg3\n")
    (d / "gl.txt").write_text("s\ng1\ng2\ng3\n")
    # deliberately wrong English name for a GO id the OBO covers
    (d / "wrong.tsv").write_text("GO:0006950\tWRONG NAME\n")
    rc = main(["-i", str(d / "net.tsv"), "--genelist", str(d / "gl.txt"),
               "--tmin", "1", "--zh", str(d / "wrong.tsv"), "-o", str(d / "out")])
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


# ---- cache keying, -f conflicts, ORA p-value shortcuts ----
def test_cache_format_and_version_are_part_of_the_key(tmp_path):
    """--format and the qenrich version must invalidate a cache, not just the mtime."""
    import json

    from qenrich import __version__

    objs = {"go": pd.DataFrame({"source": ["GO:1"], "target": ["G1"]})}
    src = wfile(tmp_path, "src.tsv", "x")
    cdir = cache_dir_for(src)
    save_objects(objs, cdir, src, "eggnog")
    assert cache_fresh(cdir, src, "eggnog")
    assert not cache_fresh(cdir, src, "generic")  # a different parse -> stale
    assert cache_fresh(cdir, src)  # fmt omitted: mtime + version only
    meta = cdir / "meta.json"
    assert json.loads(meta.read_text())["version"] == __version__
    stamp = json.loads(meta.read_text())
    stamp["version"] = "0.0.0"  # written by another qenrich
    meta.write_text(json.dumps(stamp))
    assert not cache_fresh(cdir, src, "eggnog")
    stamp.pop("version")  # legacy stamp: no version key at all
    meta.write_text(json.dumps(stamp))
    assert not cache_fresh(cdir, src)


def test_forced_format_is_not_defeated_by_the_cache(tmp_path, capsys):
    """Regression: a fresh cache used to override --format silently."""
    import json

    from qenrich._cli import main

    annot = wfile(tmp_path, "a.emapper.annotations.tsv", EGGNOG)
    gl = wfile(tmp_path, "g.txt", GENELIST)
    assert main(["-i", annot, "--genelist", gl, "--tmin", "1", "-o", str(tmp_path / "o1")]) == 0
    cdir = cache_dir_for(annot)
    assert json.loads((cdir / "meta.json").read_text())["format"] == "eggnog"
    capsys.readouterr()
    rc = main(["-i", annot, "--genelist", gl, "--format", "generic", "--tmin", "1",
               "-o", str(tmp_path / "o2")])
    assert rc == 0
    out = capsys.readouterr().out
    assert "parsed and cached" in out and "using cache" not in out
    assert json.loads((cdir / "meta.json").read_text())["format"] == "generic"


def test_no_cache_neither_reads_nor_writes(tmp_path, capsys):
    from qenrich._cli import main

    annot = wfile(tmp_path, "a.tsv", "#query\tGOs\nG1\tGO:0000001\nG2\tGO:0000001\nG3\tGO:0000001\n")
    gl = wfile(tmp_path, "g.txt", "s\nG1\nG2\nG3\n")
    rc = main(["-i", annot, "--genelist", gl, "--tmin", "1", "--no-cache", "-o", str(tmp_path / "o1")])
    assert rc == 0
    assert not cache_dir_for(annot).exists()
    out = capsys.readouterr().out
    assert "--no-cache" in out and "using cache" not in out
    assert (tmp_path / "o1" / "s_enrichment.tsv").is_file()
    assert main(["-i", annot, "--genelist", gl, "--tmin", "1", "-o", str(tmp_path / "o2")]) == 0
    assert cache_dir_for(annot).is_dir()  # a normal run does cache
    capsys.readouterr()
    assert main(["-i", annot, "--genelist", gl, "--tmin", "1", "--no-cache",
                 "-o", str(tmp_path / "o3")]) == 0
    assert "using cache" not in capsys.readouterr().out


def test_feature_conflict_on_single_net_errors(tmp_path, capsys):
    """-f naming another object must fail loudly, not enrich whatever -i pointed at."""
    from qenrich._cli import main

    annot = wfile(tmp_path, "a.tsv", "#query\tGOs\tKEGG_ko\n"
                                    "G1\tGO:0000001\tK00001\nG2\tGO:0000001\tK00002\n"
                                    "G3\tGO:0000001\tK00001\n")
    gl = wfile(tmp_path, "g.txt", "s\nG1\nG2\nG3\n")
    db = tmp_path / "db"
    assert main(["parse", annot, "-o", str(db)]) == 0
    assert main(["-i", "go", "--db", str(db), "-f", "go", "--genelist", gl, "--tmin", "1",
                 "-o", str(tmp_path / "ok")]) == 0
    capsys.readouterr()
    rc = main(["-i", "go", "--db", str(db), "-f", "kegg", "--genelist", gl, "--tmin", "1",
               "-o", str(tmp_path / "bad")])
    assert rc == 1
    err = capsys.readouterr().err
    assert "conflicts" in err and "kegg" in err
    assert not (tmp_path / "bad").exists()
    netf = wfile(tmp_path, "n.tsv", NET)
    assert main(["-i", netf, "-f", "go", "--genelist", gl, "--tmin", "1",
                 "-o", str(tmp_path / "bad2")]) == 1
    assert "conflicts" in capsys.readouterr().err


def test_fisher_pvalues_match_scipy_on_a_grid():
    """Vectorized p-values must track sts.fisher_exact for every table shape.

    Agreement is ~1e-10 relative, not bit-identical: the pmf goes through betaln
    where scipy's goes through Boost. The bound is loose enough for that rounding
    and tight enough to catch an algorithmic mistake.
    """
    from qenrich._fisher import fisher_pvalues

    tables = []
    for universe in (6, 50, 500, 3000, 30000):
        for term in (1, 3, 20, 200, universe // 3, universe):
            if term > universe:
                continue
            for size in (1, 5, 50, universe // 4, universe):
                if size > universe:
                    continue
                for k in {0, 1, term // 3, term // 2, term - 1, term, min(term, size), min(term, size) - 1}:
                    # a real term/set pair always satisfies term + size - k <= universe
                    if not max(0, term + size - universe) <= k <= min(term, size):
                        continue
                    tables.append((k, term, size, universe))
    assert len(tables) > 250
    K = np.array([t[1] for t in tables])
    n = np.array([t[2] for t in tables])
    N = np.array([t[3] for t in tables])
    obs = np.array([t[0] for t in tables])
    for alt in ("greater", "less"):
        got = fisher_pvalues(K, n, N, obs, alternative=alt)
        want = np.array([sts.fisher_exact([[k, Kt - k], [s - k, Nt - Kt - s + k]], alternative=alt)[1]
                         for k, Kt, s, Nt in tables])
        assert np.allclose(got, want, rtol=1e-9, atol=1e-12), (alt, np.abs(got - want).max())


def test_ora_never_calls_scipy_fisher(tmp_path, monkeypatch):
    """run_ora computes its p-values locally; scipy's per-table call is unused."""
    import scipy.stats as sps

    net = PARSERS["net"](wfile(tmp_path, "n.tsv", NET))["net"]
    _, sets, _ = read_genelist(wfile(tmp_path, "l.txt", GENELIST))

    def boom(*a, **k):
        raise AssertionError("run_ora must not call sts.fisher_exact")

    monkeypatch.setattr(sps, "fisher_exact", boom)
    results, _, stats = run_ora(net, sets, tmin=3)  # 2 sets x 3 terms, universe 6, term size 3
    assert stats["up"]["n_terms"] == 3
    # all three terms hold query genes, so all three are tested and reported
    up = results["up"].set_index("term")
    assert set(up.index) == {"GO:0000001", "GO:0000002", "GO:0000003"}
    assert (up["pvalue"] <= 1.0).all() and (up["pvalue"] > 0).all()


def test_terms_without_a_query_gene_are_not_tested(tmp_path):
    """clusterProfiler tests only gene sets holding a query gene; so does run_ora."""
    net = pd.DataFrame({"source": ["T1", "T1", "T2", "T2"], "target": ["g1", "g2", "g3", "g4"]})
    results, _, stats = run_ora(net, {"s": ["g1"]}, tmin=1)
    assert set(results["s"]["term"]) == {"T1"}          # T2 shares no gene -> not tested
    assert stats["s"]["n_terms"] == 1 and stats["s"]["n_empty"] == 1
    # and BH runs over the tested terms only, so padj == pvalue for a single test
    assert results["s"].iloc[0]["padj"] == results["s"].iloc[0]["pvalue"]


def test_alternative_less_tests_depletion(tmp_path):
    # universe = g0..g19 (T0 keeps g12..g19 in it); T1 covers 10 of them, T2 only 2
    net = pd.DataFrame({
        "source": ["T1"] * 10 + ["T2"] * 2 + ["T0"] * 8,
        "target": [f"g{i}" for i in range(10)] + ["g10", "g11"] + [f"g{i}" for i in range(12, 20)],
    })
    # set hits the net in {g0..g4, g10}: T1 k=5, T2 k=1, set size 6
    sets = {"s": [f"g{i}" for i in range(5)] + ["g10", "outside1", "outside2"]}
    default = run_ora(net, sets, tmin=1)[0]["s"].set_index("term")
    less = run_ora(net, sets, tmin=1, alternative="less")[0]["s"].set_index("term")
    for term, k, term_size in (("T1", 5, 10), ("T2", 1, 2)):
        table = [[k, term_size - k], [6 - k, 20 - term_size - 6 + k]]
        assert abs(default.loc[term, "pvalue"] - sts.fisher_exact(table, alternative="greater")[1]) < 1e-12
        assert abs(less.loc[term, "pvalue"] - sts.fisher_exact(table, alternative="less")[1]) < 1e-12
    # T1 is enriched, so its depletion p-value is the larger of the two
    assert less.loc["T1", "pvalue"] > default.loc["T1", "pvalue"]

    from qenrich._cli import main

    nf = wfile(tmp_path, "n.tsv", NET)
    gl = wfile(tmp_path, "g.txt", GENELIST)
    assert main(["-i", nf, "--genelist", gl, "--alternative", "less", "--no-obo", "--tmin", "3",
                 "-o", str(tmp_path / "out")]) == 0
    res = pd.read_csv(tmp_path / "out" / "up_enrichment.tsv", sep="\t")
    assert len(res) == 3  # --no-obo keeps the raw three terms; no propagation to count around


def test_ora_empty_net_returns_empty_tables():
    """An annotation net with no rows must return empty results, not raise."""
    empty = pd.DataFrame({"source": pd.Series(dtype=str), "target": pd.Series(dtype=str)})
    results, es_wide, stats = run_ora(empty, {"s": ["G1", "G2"]}, tmin=1)
    assert results["s"].empty and results["s"].columns.tolist() == [
        "term", "term_size", "overlap", "genes", "pvalue", "log_or", "padj"]
    assert stats["s"] == {"n_input": 2, "n_hit": 0, "n_terms": 0, "n_pruned": 0, "n_empty": 0}
    assert es_wide.empty


def test_ora_duplicate_pairs_count_distinct_genes():
    """A net with repeated (term, gene) rows keeps set semantics: distinct genes."""
    net = pd.DataFrame({"source": ["T1", "T1", "T1", "T2"], "target": ["g1", "g1", "g2", "g1"]})
    results, _, stats = run_ora(net, {"s": ["g1"]}, tmin=1)
    row = results["s"].set_index("term")
    assert row.loc["T1", "term_size"] == 2  # g1, g2 -> g1 counted once
    assert row.loc["T1", "overlap"] == 1 and row.loc["T1", "genes"] == "g1"
    assert stats["s"]["n_terms"] == 2


def test_tied_pvalues_keep_term_order():
    """Terms with identical tables tie exactly; their row order stays alphabetical."""
    net = pd.DataFrame({
        "source": ["T2", "T2", "T1", "T1", "T3"],
        "target": ["g1", "g2", "g1", "g2", "g1"],
    })
    # T1 and T2 have the same size and overlap with {g1} -> identical p; T3 is
    # smaller, and all three rows end up with padj == 1 here, so every one ties
    results, _, _ = run_ora(net, {"s": ["g1"]}, tmin=1)
    row = results["s"].set_index("term")
    assert row.loc["T1", "pvalue"] == row.loc["T2", "pvalue"]
    assert list(results["s"]["term"]) == ["T1", "T2", "T3"]


# ---- propagated-GO net cached next to the parsed objects ----
def test_propagated_go_cache_roundtrip(tmp_path, capsys):
    """A second run reads go_propagated.tsv instead of re-propagating."""
    from qenrich._cli import main

    annot = wfile(tmp_path, "a.tsv",
                  "#query\tGOs\nG1\tGO:0000001\nG2\tGO:0000001\nG3\tGO:0000003\nG4\tGO:0000003\n")
    gl = wfile(tmp_path, "g.txt", "s\nG1\nG2\nG3\n")
    rc = main(["-i", annot, "--genelist", gl, "--tmin", "1", "-o", str(tmp_path / "r1")])
    assert rc == 0
    cdir = tmp_path / "a.tsv.qenrich"
    assert (cdir / "go_propagated.npz").is_file() and (cdir / "go_propagated.json").is_file()
    out = capsys.readouterr().out
    assert "propagated GO DAG" in out
    first = (tmp_path / "r1" / "s_enrichment.tsv").read_text()

    rc = main(["-i", annot, "--genelist", gl, "--tmin", "1", "-o", str(tmp_path / "r2")])
    assert rc == 0
    out = capsys.readouterr().out
    assert "using propagated GO net" in out
    assert (tmp_path / "r2" / "s_enrichment.tsv").read_text() == first  # same numbers


def test_propagated_go_cache_invalidated_by_obo(tmp_path, capsys):
    """A different OBO (other file or newer mtime) must re-propagate, not serve the old net."""
    from qenrich._cli import main

    annot = wfile(tmp_path, "a.tsv",
                  "#query\tGOs\nG1\tGO:0000001\nG2\tGO:0000001\nG3\tGO:0000003\nG4\tGO:0000003\n")
    gl = wfile(tmp_path, "g.txt", "s\nG1\nG2\nG3\n")
    obo1 = wfile(tmp_path, "o1.obo", OBO)
    rc = main(["-i", annot, "--genelist", gl, "--tmin", "1", "--obo", obo1,
               "-o", str(tmp_path / "r1")])
    assert rc == 0
    capsys.readouterr()

    # same OBO, untouched mtime: cache hit
    rc = main(["-i", annot, "--genelist", gl, "--tmin", "1", "--obo", obo1,
               "-o", str(tmp_path / "r2")])
    assert "using propagated GO net" in capsys.readouterr().out

    # same file, newer mtime: stale, re-propagates
    import os
    os.utime(obo1, (2_000_000_000, 2_000_000_000))
    rc = main(["-i", annot, "--genelist", gl, "--tmin", "1", "--obo", obo1,
               "-o", str(tmp_path / "r3")])
    out = capsys.readouterr().out
    assert "propagated GO DAG" in out and "using propagated" not in out

    # and the re-stashed net is served again
    rc = main(["-i", annot, "--genelist", gl, "--tmin", "1", "--obo", obo1,
               "-o", str(tmp_path / "r4")])
    assert "using propagated GO net" in capsys.readouterr().out


def test_propagated_go_cache_respects_no_cache(tmp_path, capsys):
    """--no-cache neither reads nor writes go_propagated.tsv."""
    from qenrich._cli import main

    annot = wfile(tmp_path, "a.tsv", "#query\tGOs\nG1\tGO:0000001\nG2\tGO:0000001\n")
    gl = wfile(tmp_path, "g.txt", "s\nG1\nG2\n")
    rc = main(["-i", annot, "--genelist", gl, "--tmin", "1", "--no-cache",
               "-o", str(tmp_path / "r")])
    assert rc == 0
    out = capsys.readouterr().out
    assert "propagated GO DAG" in out and "using propagated" not in out
    assert not (tmp_path / "a.tsv.qenrich" / "go_propagated.npz").is_file()
    assert not (tmp_path / "a.tsv.qenrich").exists()  # --no-cache leaves no cache dir at all


def test_net_code_roundtrip(tmp_path):
    """save_net/load_net keep the pairs, including duplicates and empty nets."""
    from qenrich._io import load_net, save_net

    df = pd.DataFrame({"source": ["GO:1", "GO:1", "GO:2"], "target": ["g1", "g1", "g2"]})
    p = tmp_path / "n.npz"
    save_net(df, p)
    back = load_net(p)
    assert list(back.columns) == ["source", "target"]
    assert back.values.tolist() == df.values.tolist()
    assert str(back["source"].dtype) == "string"

    empty = pd.DataFrame({"source": pd.Series([], dtype="string"),
                          "target": pd.Series([], dtype="string")})
    save_net(empty, tmp_path / "e.npz")
    out = load_net(tmp_path / "e.npz")
    assert len(out) == 0 and list(out.columns) == ["source", "target"]


def test_cache_without_object_counts_is_migrated(tmp_path, capsys):
    """A stamp written before per-object counts existed is filled in, not crashed on."""
    import json
    from qenrich._cli import main

    annot = wfile(tmp_path, "a.tsv", "#query\tGOs\nG1\tGO:0000001\nG2\tGO:0000001\n")
    gl = wfile(tmp_path, "g.txt", "s\nG1\nG2\n")
    rc = main(["-i", annot, "--genelist", gl, "--tmin", "1", "-o", str(tmp_path / "r1")])
    assert rc == 0
    capsys.readouterr()
    meta_path = tmp_path / "a.tsv.qenrich" / "meta.json"
    meta = json.loads(meta_path.read_text())
    del meta["objects"]  # as an older qenrich would have stamped it
    meta_path.write_text(json.dumps(meta))

    rc = main(["-i", annot, "--genelist", gl, "--tmin", "1", "-o", str(tmp_path / "r2")])
    assert rc == 0
    out = capsys.readouterr().out
    assert "object 'go':" in out  # the listing still reports the go object
    repaired = json.loads(meta_path.read_text())["objects"]  # stamp repaired in place
    from qenrich._io import load_object
    df = load_object(tmp_path / "a.tsv.qenrich", "go")
    assert repaired["go"] == {"terms": int(df["source"].nunique()),
                              "genes": int(df["target"].nunique())}
