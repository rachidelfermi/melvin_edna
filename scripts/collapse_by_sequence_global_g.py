#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, re, sys
from typing import Dict, List, Set, Optional
import pandas as pd

# =========================================================
# Helpers
# =========================================================

# Drop any Homo (genus) and H. sapiens text patterns
HUMAN_PATTERNS = [
    r"\bHomo\s*sapiens\b",
    r"\bHomosapiens\b",
    r"s[_-]?Homo\s*sapiens",
    r"\bg[_-]?Homo\b",
    r"\bHomo\b",
]

def is_human_text(text: str) -> bool:
    if not isinstance(text, str) or not text:
        return False
    for p in HUMAN_PATTERNS:
        if re.search(p, text, flags=re.I):
            return True
    return False

def trim_ncbi_suffix(word: str) -> str:
    """Trim trailing _<digits> if present (e.g., Gallus_9030 → Gallus)."""
    return re.sub(r"_(\d+)$", "", str(word))

def to_genus_from_token(token: str) -> Optional[str]:
    """
    Normalize a genus token like 'Gallus' or 'Gallus_9030' (or 'g__Gallus') -> 'Gallus'.
    Require Genus-style capitalization (Aaaa...).
    """
    if not token:
        return None
    token = str(token)
    # remove common leading prefixes like g_, g:, g__
    token = re.sub(r"^\s*g\s*[:_]\s*", "", token, flags=re.I)
    token = token.replace("__", "_")
    token = re.sub(r"_+", " ", token).strip()
    if not token:
        return None
    # take first piece (before spaces, if any)
    g = trim_ncbi_suffix(token.split()[0])
    if re.match(r"^[A-Z][a-zA-Z-]+$", g):
        return g
    return None

# ---- Strict GENUS extraction from lineage (REQUIRED) ----
def extract_genus_from_lin_strict(lin: str) -> Optional[str]:
    """
    Extract genus strictly from lineage like '...,g_Gallus,...'.
    Requires presence of a g_ (or g:) token in the lineage.
    """
    if not isinstance(lin, str) or not lin:
        return None
    m = re.search(r"(?:\bg\s*[:_])\s*([A-Za-z][A-Za-z0-9_.-]*)", lin)
    if not m:
        return None
    return to_genus_from_token(m.group(1))

# ---- Optional loose genus extraction from tax (for cross-check only) ----
def extract_genus_from_tax_loose(tax: str) -> Optional[str]:
    """Optional: genus from SINTAX 'tax' (g:Genus...), used only to check disagreement with lin."""
    if not isinstance(tax, str) or not tax:
        return None
    m = re.search(r"(?:\bg\s*[:_])\s*([A-Za-z][A-Za-z0-9_.-]*)", tax)
    if not m:
        return None
    return to_genus_from_token(m.group(1))

def load_taxonomy(tax_tsv: str) -> pd.DataFrame:
    """
    Taxa TSV has NO HEADER:
      col0: OTU_ID
      col1: tax (SINTAX-like)
      col2+: optional; may contain lineage like 'k_/p_/.../g_/s_'
    Returns: OTU_ID, tax, lin
    """
    if not os.path.exists(tax_tsv) or os.path.getsize(tax_tsv) == 0:
        return pd.DataFrame({"OTU_ID": [], "tax": [], "lin": []})
    t = pd.read_csv(tax_tsv, sep="\t", header=None, dtype=str, na_filter=False)
    if t.shape[1] < 2:
        return pd.DataFrame({"OTU_ID": [], "tax": [], "lin": []})
    t.columns = ["OTU_ID", "tax"] + [f"extra{i}" for i in range(1, t.shape[1]-1)]
    lin_col = None
    for c in t.columns[2:]:
        blob = " ".join(t[c].head(100).astype(str).tolist()).lower()
        if any(tok in blob for tok in ["k_", "p_", "c_", "o_", "f_", "g_", "s_"]):
            lin_col = c
            break
    t["lin"] = t[lin_col] if lin_col else ""
    return t[["OTU_ID", "tax", "lin"]]

def load_fasta_ids(fa_path: str) -> Dict[str, str]:
    """FASTA → dict(OTU_ID -> sequence uppercase)"""
    seqs: Dict[str, str] = {}
    if not os.path.exists(fa_path) or os.path.getsize(fa_path) == 0:
        return seqs
    with open(fa_path) as fh:
        cur = None
        buf: List[str] = []
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = "".join(buf).upper()
    return seqs

def write_fasta(path: str, ids: List[str], id2seq: Dict[str, str]) -> None:
    with open(path, "w") as out:
        for oid in ids:
            seq = id2seq.get(oid)
            if not seq:
                continue
            out.write(f">{oid}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

# -------------------------------------------------------
# Control detection tailored to your column names
# -------------------------------------------------------

CONTROL_PREFIXES = ("ncdna", "ncpcr", "pc", "blank", "control")

def looks_like_locus_suffix(col_l: str) -> bool:
    # Your columns end with 'merged_mammal' or 'merged_12s'
    return col_l.endswith("merged_mammal") or col_l.endswith("merged_12s")

def is_control_column(col: str) -> bool:
    col_l = col.lower()
    return looks_like_locus_suffix(col_l) and col_l.startswith(CONTROL_PREFIXES)

def detect_controls(sample_names: List[str], regex: str, explicit: Set[str]) -> Set[str]:
    pat = re.compile(regex, re.I) if regex else None
    exp_l = {s.lower() for s in explicit}
    controls = set()
    for s in sample_names:
        s_l = s.lower()
        if s_l in exp_l:
            controls.add(s); continue
        if is_control_column(s):
            controls.add(s); continue
        if pat and pat.search(s):
            controls.add(s)
    return controls

# =========================================================
# STRICT genus-only base (used by both filters)
# =========================================================

def make_genus_only_base(df: pd.DataFrame, samples: List[str]) -> pd.DataFrame:
    """
    STRICT genus-level gate:
      - Keep ONLY rows where lineage (lin) has g_… (genus present).
      - genus_name is taken from lin.
      - If tax also encodes a genus and it disagrees with lin → drop row.
      - Drop Homo (genus).
    """
    df = df.copy()

    # genus from lin is REQUIRED
    df["genus_from_lin"] = df["lin"].apply(extract_genus_from_lin_strict)
    df = df.loc[df["genus_from_lin"].notna()].copy()

    # optional cross-check with 'tax' (if present)
    df["genus_from_tax"] = df["tax"].apply(extract_genus_from_tax_loose)
    mismatch = (
        df["genus_from_tax"].notna() &
        (df["genus_from_tax"] != df["genus_from_lin"])
    )
    if mismatch.any():
        df = df.loc[~mismatch].copy()

    # set canonical genus_name
    df["genus_name"] = df["genus_from_lin"]

    # drop Homo
    homo_mask = (df["genus_name"].str.fullmatch("Homo", case=False)) | \
                df["tax"].apply(is_human_text) | df["lin"].apply(is_human_text)
    df = df.loc[~homo_mask].copy()

    # tidy up helper cols
    df = df.drop(columns=["genus_from_lin", "genus_from_tax"])
    return df

# =========================================================
# Relative-only filter (genus level), mapped back to OTUs
# =========================================================

def apply_relative_only(base_df: pd.DataFrame, samples: List[str], rel_threshold: float) -> pd.DataFrame:
    counts = base_df[["genus_name"] + samples].groupby("genus_name", as_index=False).sum()
    sample_totals = {s: float(counts[s].sum()) for s in samples}
    keep = pd.DataFrame(False, index=counts["genus_name"], columns=samples)
    keep.index.name = "genus_name"  # ensure reset_index produces 'genus_name'
    counts_idx = counts.set_index("genus_name")
    for s in samples:
        tot = sample_totals[s]
        if tot > 0:
            keep.loc[counts_idx.index, s] = (counts_idx[s] / tot) > rel_threshold

    keep_long = (
        keep.reset_index()
            .melt(id_vars="genus_name", var_name="sample", value_name="keep_flag")
            .set_index(["genus_name", "sample"])["keep_flag"]
    )

    def zero_if_drop(row):
        g = row["genus_name"]
        for col in samples:
            if not keep_long.get((g, col), False):
                row[col] = 0
        return row

    df_rel = base_df.apply(zero_if_drop, axis=1)
    df_rel = df_rel.loc[df_rel[samples].sum(axis=1) > 0].copy()
    return df_rel

# =========================================================
# Control-only filter (genus level), mapped back to OTUs
# =========================================================

def apply_control_only(base_df: pd.DataFrame, samples: List[str], controls: Set[str]) -> pd.DataFrame:
    counts = base_df[["genus_name"] + samples].groupby("genus_name", as_index=False).sum()
    counts_idx = counts.set_index("genus_name")
    ctrl_cols = [s for s in samples if s in controls]
    non_ctrl_cols = [s for s in samples if s not in controls]

    keep = pd.DataFrame(False, index=counts_idx.index, columns=samples)
    keep.index.name = "genus_name"  # ensure reset_index makes this column

    if len(ctrl_cols) == 0:
        keep.loc[:, :] = True
    else:
        ctrl_max = counts_idx[ctrl_cols].max(axis=1)
        for s in non_ctrl_cols:
            keep[s] = counts_idx[s] > ctrl_max
        # Optional: zero values in control columns (typical)
        for s in ctrl_cols:
            keep[s] = False

    keep_long = (
        keep.reset_index()
            .melt(id_vars="genus_name", var_name="sample", value_name="keep_flag")
            .set_index(["genus_name", "sample"])["keep_flag"]
    )

    def zero_if_drop(row):
        g = row["genus_name"]
        for col in samples:
            if not keep_long.get((g, col), False):
                row[col] = 0
        return row

    df_ctrl = base_df.apply(zero_if_drop, axis=1)
    df_ctrl = df_ctrl.loc[df_ctrl[samples].sum(axis=1) > 0].copy()
    return df_ctrl

# =========================================================
# Writers
# =========================================================

def write_outputs(prefix: str, locus: str, df: pd.DataFrame, samples: List[str], out_dir: str, id2seq: Dict[str, str]):
    out_tab = os.path.join(out_dir, f"{prefix}_{locus}_otu_table.tsv")
    df_out = df[["OTU_ID"] + samples].copy().rename(columns={"OTU_ID": "#OTU ID"})
    df_out.to_csv(out_tab, sep="\t", index=False)

    out_rel = os.path.join(out_dir, f"{prefix}_{locus}_otu_table.relative.tsv")
    rel_df = df_out.copy()
    for s in samples:
        tot = float(rel_df[s].sum())
        rel_df[s] = (rel_df[s] / (tot if tot != 0 else 1.0)).astype(float)
    rel_df.to_csv(out_rel, sep="\t", index=False, float_format="%.8f")

    # Primary taxonomy file: exactly 3 columns, NO HEADER (compat with downstream)
    out_tax = os.path.join(out_dir, f"{prefix}_{locus}_taxa.tsv")
    df_tax3 = df[["OTU_ID", "tax", "lin"]].drop_duplicates()
    df_tax3.to_csv(out_tax, sep="\t", index=False, header=False)

    # Optional audit file with genus_name (WITH header) for your review
    out_tax_audit = os.path.join(out_dir, f"{prefix}_{locus}_taxa.with_genus.tsv")
    df_tax4 = df[["OTU_ID", "tax", "lin", "genus_name"]].drop_duplicates()
    df_tax4.to_csv(out_tax_audit, sep="\t", index=False, header=True)

    out_fa = os.path.join(out_dir, f"{prefix}_{locus}_otu.fa")
    if id2seq:
        write_fasta(out_fa, df["OTU_ID"].tolist(), id2seq)
    else:
        open(out_fa, "w").close()

    print(f"  table     -> {out_tab}")
    print(f"  relative  -> {out_rel}")
    print(f"  taxonomy  -> {out_tax}")
    print(f"  taxonomy* -> {out_tax_audit}  (includes genus_name)")
    print(f"  fasta     -> {out_fa}")

# =========================================================
# Driver per locus
# =========================================================

def run_both_filters(
    global_otutab: str,
    tax_tsv: str,
    otu_fasta: str,
    locus: str,
    out_dir: str,
    rel_threshold: float,
    control_regex: str,
    explicit_controls: List[str],
):
    os.makedirs(out_dir, exist_ok=True)

    # Load global table
    g = pd.read_csv(global_otutab, sep="\t", dtype=str)
    if "#OTU ID" not in g.columns and "OTU ID" in g.columns:
        g = g.rename(columns={"OTU ID": "#OTU ID"})
    if "#OTU ID" not in g.columns:
        raise SystemExit(f"Missing '#OTU ID' in {global_otutab}")

    samples = [c for c in g.columns if c != "#OTU ID"]
    if not samples:
        raise SystemExit(f"No sample columns in {global_otutab}")

    for c in samples:
        g[c] = pd.to_numeric(g[c], errors="coerce").fillna(0).astype(int)

    tax = load_taxonomy(tax_tsv)
    id2seq = load_fasta_ids(otu_fasta)

    # Join tax
    df = g.rename(columns={"#OTU ID": "OTU_ID"}).merge(tax, on="OTU_ID", how="left")

    # Controls
    explicit = {s.strip() for s in explicit_controls if s.strip()}
    controls = detect_controls(samples, control_regex, explicit)

    print(f"[INFO] {locus}: detected {len(controls)} control columns: "
          f"{', '.join(sorted(controls)) if controls else 'none'}")

    # Base: STRICT genus-only + drop Homo
    base_df = make_genus_only_base(df, samples)

    # ---- Relative-only outputs ----
    print(f"[OK] {locus} RELATIVE-ONLY (thr={rel_threshold}):")
    rel_df = apply_relative_only(base_df, samples, rel_threshold)
    write_outputs(prefix="filtered_rel", locus=locus, df=rel_df, samples=samples, out_dir=out_dir, id2seq=id2seq)

    # ---- Control-only outputs ----
    print(f"[OK] {locus} CONTROL-ONLY:")
    ctrl_df = apply_control_only(base_df, samples, controls)
    write_outputs(prefix="filtered_control", locus=locus, df=ctrl_df, samples=samples, out_dir=out_dir, id2seq=id2seq)

# =========================================================
# CLI
# =========================================================

def main():
    p = argparse.ArgumentParser(
        description=(
            "STRICT genus-only (from LIN), drop Homo, then produce TWO output sets per locus:\n"
            "  • filtered_rel_*     = relative-only at genus level (thr=--rel-threshold)\n"
            "  • filtered_control_* = control-only at genus level (sample count > max control)\n"
            "Outputs mapped back to OTU IDs; no collapsing."
        )
    )
    p.add_argument("--input-dir",  default=".", help="Dir with *_otu_table.tsv, *_taxa.tsv, *_otu.fa")
    p.add_argument("--out-dir",    default=".", help="Where to write filtered outputs")
    p.add_argument("--rel-threshold", type=float, default=0.01,
                   help="Per-sample relative abundance cutoff at GENUS level (default 0.01 = 1%)")
    p.add_argument("--control-regex", default=r"",
                   help="Extra regex to mark control columns (case-insensitive); applied in addition to built-in prefixes")
    p.add_argument("--control-list", default="",
                   help="Comma-separated explicit control column names (exact matches)")
    args = p.parse_args()

    IN  = args.input_dir
    OUT = args.out_dir
    rel_thr = float(args.rel_threshold)
    ctrl_regex = args.control_regex
    ctrl_list = [x.strip() for x in args.control_list.split(",")] if args.control_list else []

    # inputs per locus
    f12_tab, f12_tax, f12_fa = os.path.join(IN, "12S_otu_table.tsv"),     os.path.join(IN, "12S_taxa.tsv"),     os.path.join(IN, "12S_otu.fa")
    fm_tab,  fm_tax,  fm_fa  = os.path.join(IN, "mammal_otu_table.tsv"),  os.path.join(IN, "mammal_taxa.tsv"),  os.path.join(IN, "mammal_otu.fa")

    if os.path.exists(f12_tab) and os.path.exists(f12_tax) and os.path.exists(f12_fa):
        run_both_filters(f12_tab, f12_tax, f12_fa, "12S", OUT, rel_thr, ctrl_regex, ctrl_list)
    else:
        print("[WARN] Skipping 12S (missing one or more inputs)", file=sys.stderr)

    if os.path.exists(fm_tab) and os.path.exists(fm_tax) and os.path.exists(fm_fa):
        run_both_filters(fm_tab, fm_tax, fm_fa, "mammal", OUT, rel_thr, ctrl_regex, ctrl_list)
    else:
        print("[WARN] Skipping mammal (missing one or more inputs)", file=sys.stderr)

if __name__ == "__main__":
    main()
