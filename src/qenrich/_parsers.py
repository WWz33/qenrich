"""Parsers: each mainstream annotation format to ``{object_name: net DataFrame}``.

Every net DataFrame has columns ``source`` (term/family) and ``target`` (gene),
optionally ``weight``. ORA ignores weights; they are kept only when the source
format provides a meaningful one (none does for v1).
"""

import re
import warnings

import pandas as pd

from ._io import open_text
from ._sniff import GO_RE, IPR_RE, KO_RE, PFAM_RE

_SKIP = {"", "-", "NA", "N/A", "None"}


def _read_tsv(path: str, **kw) -> pd.DataFrame:
    """Read a TSV; ragged rows are truncated, not turned into an index column
    (pandas' default would silently shift every column by one)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", pd.errors.ParserWarning)
        return pd.read_csv(path, sep="\t", index_col=False, **kw)


def _pairs(rows: list[tuple[str, str]]) -> pd.DataFrame:
    df = pd.DataFrame(sorted(set(rows)), columns=["source", "target"])
    return df[~df["source"].isin(_SKIP) & ~df["target"].isin(_SKIP)]


def _split_ids(s: str) -> list[str]:
    return [t for t in re.split(r"[|,;\s]+", str(s).strip()) if t and t not in _SKIP]


def _net(df: pd.DataFrame, source: str, target: str) -> pd.DataFrame:
    out = df[[source, target]].copy()
    out.columns = ["source", "target"]
    return out.drop_duplicates().reset_index(drop=True)


def _header_index(path, pred) -> int:
    """Line index of first line satisfying pred (used to locate header rows)."""
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            if pred(line):
                return i
    raise ValueError(f"header row not found in {path}")


# ---------------------------------------------------------------- eggNOG ----
_COG = {
    "A": "RNA processing and modification", "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion", "D": "Cell cycle control and division",
    "E": "Amino acid transport and metabolism", "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism", "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism", "J": "Translation and ribosomal biogenesis",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane biogenesis", "N": "Cell motility",
    "O": "Posttranslational modification and chaperones", "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis", "R": "General function prediction only",
    "S": "Function unknown", "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking and secretion", "V": "Defense mechanisms",
    "W": "Extracellular structures", "Y": "Nuclear structure", "Z": "Cytoskeleton",
}


def parse_eggnog(path: str, annot_lvl: str | int | None = None) -> dict[str, pd.DataFrame]:
    h = _header_index(path, lambda l: l.lower().startswith("#query\t") or l.lower().startswith("query\t"))
    # index_col=False: a trailing-tab row would otherwise make pandas eat 'query'
    # as the index and silently shift every column by one
    df = _read_tsv(path, skiprows=h, dtype=str, keep_default_na=False)
    df.columns = [c.lstrip("#").strip().lower() for c in df.columns]  # header case varies
    lvl_col = next((c for c in ("max_annot_lvl", "tax_ceiling") if c in df.columns), None)
    if annot_lvl is not None and lvl_col:
        # emapper>=2.1.13 stores "NNN|Taxon"; match the taxid before the separator
        lvl = df[lvl_col].astype(str)
        df = df[(lvl == str(annot_lvl)) | (lvl.str.split("|").str[0] == str(annot_lvl))]
    gene = df["query"].astype(str)
    objects: dict[str, pd.DataFrame] = {}

    def explode(col: str, tfx) -> None:
        if col not in df.columns:
            return
        rows = []
        for g, val in zip(gene, df[col], strict=True):
            for t in tfx(val):
                rows.append((t, g))
        if rows:
            objects[_objname(col)] = _pairs(rows)

    explode("gos", lambda v: GO_RE.findall(v))
    explode("kegg_ko", lambda v: [t for t in re.findall(r"K\d{5}", v)])
    explode("kegg_pathway", lambda v: _split_ids(v))
    # pfams holds domain *names* ("Profilin") and sometimes PF accessions; keep both
    explode("pfams", lambda v: [re.sub(r"^(PF\d{5})\.\d+$", r"\1", t) for t in _split_ids(v)])
    explode("cog_category", lambda v: [(f"{c} {_COG.get(c, '')}").strip() for c in str(v) if c in _COG])
    return objects


def _objname(col: str) -> str:
    return {"gos": "go", "kegg_ko": "kegg", "kegg_pathway": "pathway", "pfams": "pfam", "cog_category": "cog"}[col]


# ----------------------------------------------------------- InterProScan ----
def parse_iprscan(path: str) -> dict[str, pd.DataFrame]:
    df = _read_tsv(path, header=None, dtype=str, keep_default_na=False)
    ncols = df.shape[1]
    gene = df[0]
    objects: dict[str, pd.DataFrame] = {}

    def build(col: int, tfx, mask=None) -> list[tuple[str, str]]:
        rows = []
        sub = df if mask is None else df[mask]
        for g, val in zip(gene[sub.index], sub[col], strict=True):
            for t in tfx(val):
                rows.append((t, g))
        return rows

    if ncols >= 14:  # GO column (with --goterms)
        rows = build(13, lambda v: [GO_RE.search(t).group() for t in v.split(",") if GO_RE.search(t)])
        if rows:
            objects["go"] = _pairs(rows)
    if ncols >= 12:  # InterPro accession
        rows = build(11, lambda v: [t for t in v.split("|") if IPR_RE.fullmatch(t.strip())])
        if rows:
            objects["interpro"] = _pairs(rows)
    rows = build(4, lambda v: [re.sub(r"\.\d+$", "", t) for t in re.findall(r"PF\d{5}(?:\.\d+)?", v)], mask=df[3] == "Pfam")
    if rows:
        objects["pfam"] = _pairs(rows)
    return objects


# -------------------------------------------------------------- Blast2GO ----
def parse_b2g_annot(path: str) -> dict[str, pd.DataFrame]:
    go_rows, ipr_rows = [], []
    with open_text(path) as fh:
        for line in fh:
            r = [c for c in line.rstrip("\r\n\t ").split("\t")]
            if len(r) < 4 or r[0].startswith("!"):
                continue
            go_rows += [(t, r[0]) for t in GO_RE.findall(r[1])]
            # column 4 is the InterPro accession; extra columns (descriptions) are ignored
            ipr_rows += [(t, r[0]) for t in IPR_RE.findall(r[3].strip())]
    objects = {}
    if go_rows:
        objects["go"] = _pairs(go_rows)
    if ipr_rows:
        objects["interpro"] = _pairs(ipr_rows)
    return objects


def parse_b2g_tabular(path: str) -> dict[str, pd.DataFrame]:
    df = _read_tsv(path, dtype=str, keep_default_na=False)
    low = {c: c.strip().lower() for c in df.columns}
    gene_col = next((c for c, l in low.items() if l == "sequence name"), df.columns[0])
    go_col = next((c for c in df.columns if c != gene_col and df[c].astype(str).str.contains(GO_RE).any()), None)
    if go_col is None:
        raise ValueError(f"no GO column found in {path}")
    rows = [(m.group(), g) for g, v in zip(df[gene_col], df[go_col]) for m in GO_RE.finditer(str(v))]
    return {"go": _pairs(rows)} if rows else {}


# --------------------------------------------------------------- PANNZER ----
def parse_pannzer(path: str) -> dict[str, pd.DataFrame]:
    df = _read_tsv(path, dtype=str, keep_default_na=False)
    df.columns = [c.strip().lstrip("#").lower() for c in df.columns]
    gene_col = "qseqid" if "qseqid" in df.columns else "qpid"
    go_col = "go_id" if "go_id" in df.columns else next(
        (c for c in df.columns if df[c].astype(str).str.contains(GO_RE).any()), None)
    if go_col is None:
        raise ValueError(f"no GO column found in {path}")
    rows = [(m.group(), g) for g, v in zip(df[gene_col], df[go_col]) for m in GO_RE.finditer(str(v))]
    return {"go": _pairs(rows)} if rows else {}


# ------------------------------------------------------------- Trinotate ----
def parse_trinotate(path: str) -> dict[str, pd.DataFrame]:
    df = _read_tsv(path, dtype=str, keep_default_na=False)
    df.columns = [c.strip().lstrip("#").lower() for c in df.columns]
    gene_col = "gene_id" if "gene_id" in df.columns else df.columns[0]
    go_cols = [c for c in df.columns if "gene_ontology" in c]
    if not go_cols:
        raise ValueError(f"no gene_ontology_* columns in {path}")
    rows = []
    for _, r in df.iterrows():
        for c in go_cols:
            for m in GO_RE.finditer(str(r[c])):  # pipe- and comma-lists both allowed
                rows.append((m.group(), str(r[gene_col])))
    return {"go": _pairs(rows)} if rows else {}


# ------------------------------------------------------------------- GAF ----
def parse_gaf(path: str) -> dict[str, pd.DataFrame]:
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            r = line.rstrip("\r\n").split("\t")
            if len(r) < 15:
                continue
            if r[3].startswith("NOT"):
                continue
            rows += [(m, r[1]) for m in GO_RE.findall(r[4])]  # col 5 may be pipe-separated
    if not rows:
        raise ValueError(f"no usable annotations in GAF file: {path}")
    return {"go": _pairs(rows)}


# ------------------------------------------------------------------ GFF3 ----
def parse_gff3(path: str) -> dict[str, pd.DataFrame]:
    entries: list[tuple[str, list[str], list[str]]] = []
    parents: dict[str, str] = {}

    def attrs(s: str) -> dict[str, str]:
        out = {}
        for kv in s.split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                out[k.strip()] = v.strip()
        return out

    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            r = [c for c in line.rstrip("\r\n\t ").split("\t")]
            if len(r) != 9:
                continue
            a = attrs(r[8])
            gid = a.get("ID")
            if not gid:
                continue
            if a.get("Parent"):
                parents[gid] = a["Parent"].split(",")[0]  # first parent for fused loci
            go = GO_RE.findall(a.get("Ontology_term", ""))
            ipr = re.findall(r"IPR\d{6}", a.get("Dbxref", ""))
            if go or ipr:
                entries.append((gid, go, ipr))

    def top_id(gid: str) -> str:
        # Roll child-feature annotations (mRNA etc.) up to the top-level feature,
        # then drop namespace prefixes like "gene:" so IDs match gene lists.
        seen: set[str] = set()
        while gid in parents and gid not in seen:
            seen.add(gid)
            gid = parents[gid]
        return gid.split(":")[-1]

    go_rows, ipr_rows = [], []
    for gid, go, ipr in entries:
        t = top_id(gid)
        go_rows += [(term, t) for term in go]
        ipr_rows += [(term, t) for term in ipr]
    objects = {}
    if go_rows:
        objects["go"] = _pairs(go_rows)
    if ipr_rows:
        objects["interpro"] = _pairs(ipr_rows)
    return objects


# ------------------------------------------------------------ KofamKOALA ----
def parse_kofam(path: str) -> dict[str, pd.DataFrame]:
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            r = line.rstrip("\r\n").split("\t")
            if len(r) < 2 or r[1] == "-":
                continue
            ko = r[1].strip()
            if KO_RE.fullmatch(ko):  # reject pipe-joined cells, validate one KO per row
                rows.append((ko, r[0].strip()))
            else:
                for k in _split_ids(r[1]):  # tolerate separator-joined KOs
                    if KO_RE.fullmatch(k):
                        rows.append((k, r[0].strip()))
    if not rows:
        raise ValueError(f"no KO assignments in {path}")
    return {"kegg": _pairs(rows)}


# --------------------------------------------------------------- generic ----
_PATTERNS = [("go", GO_RE), ("kegg", KO_RE), ("pfam", PFAM_RE), ("interpro", IPR_RE)]


def parse_generic(path: str) -> dict[str, pd.DataFrame]:
    rows_raw = [line.rstrip("\r\n").split("\t") if "\t" in line else line.split() for line in open_text(path) if line.strip()]
    ncols = max(len(r) for r in rows_raw)

    def idlike(cells: list[str]) -> list[str]:
        return [t for c in cells for t in _split_ids(c)]

    first_vals = idlike(rows_raw[0][1:ncols])
    # row 0 carrying ID-like values means there is no header; dropping it here
    # silently lost the first gene of header-less files
    data = rows_raw[1:] if (not first_vals and len(rows_raw) >= 2) else rows_raw
    objects: dict[str, list[tuple[str, str]]] = {n: [] for n, _ in _PATTERNS}
    for r in data:
        r = r + [""] * (ncols - len(r))
        for j in range(1, ncols):
            for name, rex in _PATTERNS:
                for t in _split_ids(r[j]):
                    if rex.fullmatch(t):
                        objects[name].append((t, r[0].strip()))
    out = {n: _pairs(v) for n, v in objects.items() if v}
    if not out:
        raise ValueError(f"no GO/KO/PFAM/IPR identifiers found in {path}")
    return out


# ------------------------------------------------------------------- net ----
def parse_net(path: str) -> dict[str, pd.DataFrame]:
    df = _read_tsv(path, dtype=str, keep_default_na=False)
    df.columns = [c.strip().lower() for c in df.columns]
    if "source" not in df.columns or "target" not in df.columns:
        raise ValueError(f"net file must have 'source' and 'target' columns: {path}")
    return {"net": _net(df, "source", "target")}


PARSERS = {
    "eggnog": parse_eggnog,
    "iprscan": parse_iprscan,
    "b2g_annot": parse_b2g_annot,
    "b2g_tabular": parse_b2g_tabular,
    "pannzer": parse_pannzer,
    "trinotate": parse_trinotate,
    "gaf": parse_gaf,
    "gff3": parse_gff3,
    "kofam": parse_kofam,
    "generic": parse_generic,
    "net": parse_net,
}
