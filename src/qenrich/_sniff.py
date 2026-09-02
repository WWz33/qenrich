"""Format sniffing for annotation files, by fixed priority order."""

import re
from pathlib import Path

from ._io import open_text

GO_RE = re.compile(r"GO:\d{7}")
KO_RE = re.compile(r"^K\d{5}$")
PFAM_RE = re.compile(r"^PF\d{5}$")
IPR_RE = re.compile(r"^IPR\d{6}$")
MD5_RE = re.compile(r"^[0-9a-f]{32}$")

FORMATS = [
    "net",
    "eggnog",
    "gff3",
    "iprscan",
    "gaf",
    "b2g_tabular",
    "trinotate",
    "pannzer",
    "b2g_annot",
    "kofam",
    "generic",
]


def _rows(path: str | Path, n: int = 100) -> list[list[str]]:
    """First n non-empty lines, tab-split (whitespace-split if no tabs)."""
    out = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            out.append(line.rstrip("\r\n").split("\t") if "\t" in line else line.split())
            if len(out) >= n:
                break
    return out


def _cols(rows: list[list[str]]) -> int:
    return max(len(r) for r in rows)


def _header_cells(rows: list[list[str]]) -> list[list[str]]:
    """Header-candidate rows: leading-# rows plus the first body row.

    eggNOG puts its real header after ``##`` comment lines, so scanning only
    row 0 misses it.
    """
    out = []
    for r in rows:
        if r[0].startswith("#"):
            out.append(r)
        else:
            out.append(r)
            break
    return out


def _has_cols(rows: list[list[str]], names: set[str]) -> bool:
    for r in _header_cells(rows):
        cells = {c.strip().lstrip("#").lower() for c in r}
        if names <= cells:
            return True
    return False


def _is_net(rows: list[list[str]]) -> bool:
    return _has_cols(rows, {"source", "target"})


def _is_eggnog(rows: list[list[str]]) -> bool:
    return _has_cols(rows, {"query", "gos"})


def _is_gff3(rows: list[list[str]]) -> bool:
    if rows and rows[0] and rows[0][0].startswith("##gff-version"):
        return True
    return any(len(r) == 9 and any("Ontology_term=" in c for c in r[8:]) for r in rows if len(r) >= 9)


def _is_iprscan(rows: list[list[str]]) -> bool:
    body = [r for r in rows if not r[0].startswith("#")]
    if len(body) < 2 or _cols(body) < 12:
        return False
    return all(MD5_RE.match(r[1]) and r[2].isdigit() for r in body if len(r) >= 3)


def _is_gaf(rows: list[list[str]]) -> bool:
    body = [r for r in rows if not r[0].startswith("!")]
    if not body:
        return False
    wide = [r for r in body if len(r) >= 15]
    if not wide:
        return False
    hits = sum(1 for r in wide if GO_RE.search(r[4]) and re.match(r"^[A-Z]{2,4}$", r[6]))
    return hits >= max(1, len(wide) // 2)


def _is_b2g_tabular(rows: list[list[str]]) -> bool:
    return any(
        "sequence name" in {c.strip().lower() for c in r} and any("go" in c.strip().lower() for c in r)
        for r in _header_cells(rows)
    )


def _is_trinotate(rows: list[list[str]]) -> bool:
    return _has_cols(rows, {"gene_id", "transcript_id"})


def _is_pannzer(rows: list[list[str]]) -> bool:
    return _has_cols(rows, {"qpid"}) or _has_cols(rows, {"qseqid"})


def _is_b2g_annot(rows: list[list[str]]) -> bool:
    body = [r for r in rows if not r[0].startswith("!")]
    four = [r for r in body if len(r) == 4]
    if len(four) < max(2, len(body) // 2):
        return False
    return all(r[2].lower() == "interpro" for r in four)


def _is_kofam(rows: list[list[str]]) -> bool:
    body = [r for r in rows if not r[0].startswith("#")]
    if len(body) < 2 or any(len(r) < 2 for r in body):
        return False
    return all(KO_RE.match(r[1]) or r[1] == "-" for r in body)


def _is_generic(rows: list[list[str]]) -> bool:
    return _cols(rows) >= 2


_CHECKS = [
    ("net", _is_net),
    ("eggnog", _is_eggnog),
    ("gff3", _is_gff3),
    ("iprscan", _is_iprscan),
    ("gaf", _is_gaf),
    ("b2g_tabular", _is_b2g_tabular),
    ("trinotate", _is_trinotate),
    ("pannzer", _is_pannzer),
    ("b2g_annot", _is_b2g_annot),
    ("kofam", _is_kofam),
    ("generic", _is_generic),
]


def sniff(path: str | Path) -> str:
    """Return detected format name; ``generic`` is the guaranteed fallback."""
    rows = _rows(path)
    if not rows:
        raise ValueError(f"empty file: {path}")
    for fmt, check in _CHECKS:
        if check(rows):
            return fmt
    raise ValueError(f"cannot recognize format of {path}; use --format to force one of {FORMATS}")


# Human-readable names for CLI output
FORMAT_LABELS = {
    "net": "standard net TSV (source/target)",
    "eggnog": "eggNOG-mapper annotations",
    "gff3": "GFF3 with Ontology_term attributes",
    "iprscan": "InterProScan TSV",
    "gaf": "GAF 2.x gene association",
    "b2g_tabular": "Blast2GO/OmicsBox tabular export",
    "trinotate": "Trinotate report",
    "pannzer": "PANNZER2 output",
    "b2g_annot": "Blast2GO .annot",
    "kofam": "KofamKOALA/GhostKOALA/KAAS KO annotation",
    "generic": "generic table (gene + ID columns)",
}
