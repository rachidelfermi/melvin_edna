#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report.py — per-locus sample report with consistent 99% mapping.

Outputs (in --input-dir):
  sample_progress_mammal.tsv
  sample_progress_12S.tsv

Columns:
  Sample, MergedReads, FASTAseqs, DerepSeqs, UNOISE_ASVs, Mapped99, ASVs, ASVs_withTaxa
"""

import os, sys, csv, argparse, glob
from typing import Dict, List, Tuple, Set

# ---------- utils ----------
def newest_match(pats: List[str]) -> str:
    cand=[]
    for pat in pats:
        for p in glob.glob(pat):
            try: m=os.path.getmtime(p)
            except: m=0
            cand.append((m,p))
    if not cand: return ""
    cand.sort(reverse=True)
    return cand[0][1]

def find_reports_dir(results_dir: str) -> str:
    parent = os.path.abspath(os.path.join(results_dir, os.pardir))
    rep = newest_match([
        os.path.join(parent, "07_reports*"),
        os.path.join(parent, "*07_reports*")
    ])
    if not rep:
        sys.stderr.write(f"❌ Could not locate 07_reports* next to {results_dir}\n")
        sys.exit(2)
    return rep

def read_simple_counts(path: str, key_name="Sample", val_idx=1) -> Dict[str,int]:
    d={}
    if not (path and os.path.exists(path)): return d
    with open(path, newline="") as f:
        r=csv.reader(f, delimiter="\t"); header=next(r, None)
        for row in r:
            if not row or row[0] in ("Total","",key_name): continue
            try: d[row[0]]=int(float(row[val_idx]))
            except: d[row[0]]=0
    return d

def read_step9(path: str) -> Dict[str,Dict[str,int]]:
    D={}
    if not (path and os.path.exists(path)): return D
    with open(path, newline="") as f:
        r=csv.reader(f, delimiter="\t"); header=next(r, None)
        for row in r:
            if not row or row[0]=="Total": continue
            smp=row[0]
            vals=[]
            for i in range(1,4):
                try: vals.append(int(float(row[i])))
                except: vals.append(0)
            D[smp]={"FASTAseqs":vals[0], "DerepSeqs":vals[1], "UNOISE_ASVs":vals[2]}
    return D

def read_mapped_sums(otutab_path: str) -> Dict[str,int]:
    """Sum of mapped reads per sample from the 99% OTU table."""
    sums={}
    if not (otutab_path and os.path.exists(otutab_path)): return sums
    with open(otutab_path, newline="") as f:
        r=csv.reader(f, delimiter="\t"); header=next(r, None)
        if not header: return sums
        samples=header[1:]; sums={s:0 for s in samples}
        for row in r:
            if not row: continue
            for i,s in enumerate(samples, start=1):
                try: sums[s]+=int(float(row[i]))
                except: pass
    return sums

def compute_asv_presence_counts(otutab_path: str, taxa_path: str) -> Tuple[Dict[str,int], Dict[str,int]]:
    """
    From the 99% OTU table:
      - ASVs per sample = # OTUs with count>0 in that sample
      - ASVs_withTaxa per sample = those OTUs intersected with OTUs that have non-empty lineage in taxa.tsv
    Assumes taxa.tsv first column is the remapped OTU ID (OTU_*), with lineage in column 4.
    """
    if not (otutab_path and os.path.exists(otutab_path)):
        return {}, {}
    assigned: Set[str] = set()
    if taxa_path and os.path.exists(taxa_path):
        with open(taxa_path, newline="") as f:
            r=csv.reader(f, delimiter="\t")
            for row in r:
                if not row: continue
                otu=row[0]
                lineage=row[3] if len(row)>3 else ""
                if otu and lineage and lineage.strip():
                    assigned.add(otu)

    asv_counts: Dict[str,int]={}
    tax_counts: Dict[str,int]={}
    with open(otutab_path, newline="") as f:
        r=csv.reader(f, delimiter="\t")
        header=next(r, None)
        if not header: return {}, {}
        samples=header[1:]
        # track presence per sample
        present = {s:set() for s in samples}
        for row in r:
            if not row: continue
            otu=row[0]
            for i,s in enumerate(samples, start=1):
                try:
                    if int(float(row[i]))>0:
                        present[s].add(otu)
                except:
                    pass
        for s in samples:
            asv_counts[s]=len(present[s])
            tax_counts[s]=sum(1 for otu in present[s] if otu in assigned)
    return asv_counts, tax_counts

def write_tsv(path: str, header: List[str], rows: List[List]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path,"w",newline="") as f:
        w=csv.writer(f, delimiter="\t"); w.writerow(header); w.writerows(rows)

# ---------- per-locus report ----------
def build_report(locus: str, reports_dir: str, results_dir: str) -> str:
    if locus=="mammal":
        merged_path = os.path.join(reports_dir, "step6_merged_counts_mammal.tsv")
        step9_path  = os.path.join(reports_dir, "step9_per_sample_asv_mammal.tsv")
        otutab_path = os.path.join(results_dir, "mammal_otu_table.tsv")
        taxa_path   = os.path.join(results_dir, "mammal_taxa.tsv")
        out_path    = os.path.join(results_dir, "sample_progress_mammal.tsv")
    else:
        merged_path = os.path.join(reports_dir, "step6_merged_counts_12S.tsv")
        step9_path  = os.path.join(reports_dir, "step9_per_sample_asv_12S.tsv")
        otutab_path = os.path.join(results_dir, "12S_otu_table.tsv")
        taxa_path   = os.path.join(results_dir, "12S_taxa.tsv")
        out_path    = os.path.join(results_dir, "sample_progress_12S.tsv")

    merged = read_simple_counts(merged_path)
    step9  = read_step9(step9_path)
    mapped = read_mapped_sums(otutab_path)              # 99% sums
    asv_counts, tax_counts = compute_asv_presence_counts(otutab_path, taxa_path)  # from same 99% table

    samples = sorted(set(merged)|set(step9)|set(mapped)|set(asv_counts)|set(tax_counts))
    header=["Sample","MergedReads","FASTAseqs","DerepSeqs","UNOISE_ASVs","Mapped99","ASVs","ASVs_withTaxa"]
    rows=[]
    for s in samples:
        rows.append([
            s,
            merged.get(s,0),
            step9.get(s,{}).get("FASTAseqs",0),
            step9.get(s,{}).get("DerepSeqs",0),
            step9.get(s,{}).get("UNOISE_ASVs",0),
            mapped.get(s,0),
            asv_counts.get(s,0),
            tax_counts.get(s,0),
        ])
    write_tsv(out_path, header, rows)
    return out_path

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Per-locus report with consistent 99% mapping.")
    ap.add_argument("--input-dir", required=True, help="Path to 14_results* directory")
    args = ap.parse_args()

    RESULTS_DIR = os.path.abspath(args.input_dir)
    if not os.path.isdir(RESULTS_DIR):
        sys.stderr.write(f"❌ --input-dir not found: {RESULTS_DIR}\n"); sys.exit(2)
    REPORTS_DIR = find_reports_dir(RESULTS_DIR)

    print(f"ℹ️ RESULTS_DIR = {RESULTS_DIR}")
    print(f"ℹ️ REPORTS_DIR = {REPORTS_DIR}")

    m_path = build_report("mammal", REPORTS_DIR, RESULTS_DIR)
    s_path = build_report("12S", REPORTS_DIR, RESULTS_DIR)

    print(f"✅ Wrote: {m_path}")
    print(f"✅ Wrote: {s_path}")
    # NOTE: deliberately NOT writing a combined file anymore.

    return 0

if __name__ == "__main__":
    sys.exit(main())
