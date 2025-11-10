#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --account=omics_core-rqyo8fdrbpw-LOW-CPU
#SBATCH --qos=low-cpu
#SBATCH --job-name=melvin_from_06merged_AtoZ_no_minimap

set -euo pipefail
shopt -s nullglob dotglob

# === Env ===
eval "$(/srv/data/rachid.elfermi/miniconda3/bin/conda shell.bash hook)"
conda activate bioinfo

# =========================
# CONFIG
# =========================
# Prefer user-specified THREADS; otherwise fall back to SLURM, then 12
if [[ -z "${THREADS:-}" ]]; then
  if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    THREADS="$SLURM_CPUS_PER_TASK"
  elif [[ -n "${SLURM_NTASKS:-}" ]]; then
    THREADS="$SLURM_NTASKS"
  else
    THREADS=12
  fi
fi

# Fixed suffix (no host filtering)
SUF="_no_minimap"

WORKDIR="${WORKDIR:-/home/rachid.elfermi/lustre/omics_core-rqyo8fdrbpw/users/rachid.elfermi/dada2_nfcore/melvin}"
cd "$WORKDIR" || { echo "❌ cd fail: $WORKDIR"; exit 1; }

MERGED_DIR="$WORKDIR/06_merged"                  # input: mammal/ and 12S/

# Outputs (parallel structure with suffix)
HOSTFREE_ROOT="$WORKDIR/08_hostfree${SUF}"
PER_ASV_DIR="$WORKDIR/09_per_sample_asv${SUF}"
GLOBAL_DIR="$WORKDIR/09_global_merge${SUF}"
CHIMERA_DIR="$WORKDIR/10_chimeras${SUF}"
MAP_DIR="$WORKDIR/11_read_mapping${SUF}"
MIDORI_DIR="$WORKDIR/MIDORI_db"                  # unchanged
TAX_DIR="$WORKDIR/13_taxonomy${SUF}"
RESULTS_DIR="$WORKDIR/14_results${SUF}"
REPORTS_DIR="$WORKDIR/07_reports${SUF}"
SCRIPTS_DIR="$WORKDIR/scripts"
mkdir -p "$HOSTFREE_ROOT"/{mammal,12S} \
         "$PER_ASV_DIR"/{mammal,12S} \
         "$GLOBAL_DIR"/{mammal,12S} \
         "$CHIMERA_DIR" "$MAP_DIR"/{mammal,12S} \
         "$MIDORI_DIR" "$TAX_DIR" "$RESULTS_DIR" "$REPORTS_DIR"

# Export for downstream Python blocks
export RESULTS_DIR REPORTS_DIR MAP_DIR GLOBAL_DIR CHIMERA_DIR

# MIDORI UNIQ (SINTAX)
MIDORI_FASTA="$MIDORI_DIR/MIDORI2_UNIQ_NUC_GB267_srRNA_SINTAX.fasta"
MIDORI_UDB="$MIDORI_DIR/MIDORI2_UNIQ_NUC_GB267_srRNA_SINTAX.udb"
SINTAX_CUTOFF=0.8

# =========================
# HELPERS
# =========================
count_fastq_reads () {
  local f="$1"
  if [[ -f "$f" ]]; then
    if [[ "$f" == *.gz ]]; then
      zcat -f -- "$f" 2>/dev/null | grep -c '^+$' || echo 0
    else
      grep -c '^+$' "$f" 2>/dev/null || echo 0
    fi
  else
    echo 0
  fi
}

tsv_total_row () {
  local tsv="$1"
  [[ -s "$tsv" ]] || return 0
  awk 'BEGIN{FS=OFS="\t"} NR==1{next} {for(i=2;i<=NF;i++) S[i]+=$i}
       END{printf "Total"; for(i=2;i<=NF;i++) printf "\t%d", S[i]; printf "\n"}' "$tsv" >> "$tsv"
}

need() { [[ -s "$1" ]] || { echo "❌ Missing: $1"; exit 1; }; }

# =========================
# STEP 6 — Baseline counts on merged inputs
# =========================
echo -e "Sample\tMergedReads_mammal" > "$REPORTS_DIR/step6_merged_counts_mammal.tsv"
for fq in "$MERGED_DIR/mammal/"*merged_mammal.fastq; do
  s=$(basename "$fq" .fastq)
  echo -e "$s\t$(count_fastq_reads "$fq")" >> "$REPORTS_DIR/step6_merged_counts_mammal.tsv"
done
tsv_total_row "$REPORTS_DIR/step6_merged_counts_mammal.tsv"

echo -e "Sample\tMergedReads_12S" > "$REPORTS_DIR/step6_merged_counts_12S.tsv"
for fq in "$MERGED_DIR/12S/"*merged_12S.fastq; do
  s=$(basename "$fq" .fastq)
  echo -e "$s\t$(count_fastq_reads "$fq")" >> "$REPORTS_DIR/step6_merged_counts_12S.tsv"
done
tsv_total_row "$REPORTS_DIR/step6_merged_counts_12S.tsv"

# =========================
# STEP 6.5 — copy merged fastqs → “hostfree” filenames (no host report)
# =========================
while IFS=$'\t' read -r sample _; do
  [[ "$sample" == "Sample" || -z "$sample" || "$sample" == "Total" ]] && continue
  src="$MERGED_DIR/mammal/${sample}.fastq"
  out="$HOSTFREE_ROOT/mammal/${sample}.hostfree.fq"
  [[ -s "$src" ]] && cp -f "$src" "$out" || :
done < "$REPORTS_DIR/step6_merged_counts_mammal.tsv"

while IFS=$'\t' read -r sample _; do
  [[ "$sample" == "Sample" || -z "$sample" || "$sample" == "Total" ]] && continue
  src="$MERGED_DIR/12S/${sample}.fastq"
  out="$HOSTFREE_ROOT/12S/${sample}.hostfree.fq"
  [[ -s "$src" ]] && cp -f "$src" "$out" || :
done < "$REPORTS_DIR/step6_merged_counts_12S.tsv"

# =========================
# STEP 7–9 — PER-SAMPLE ASV CALLING (FASTA → derep → UNOISE)
# =========================
echo -e "Sample\tFASTAseqs\tDerepSeqs\tUNOISE_ASVs" > "$REPORTS_DIR/step9_per_sample_asv_mammal.tsv"
for fq in "$HOSTFREE_ROOT/mammal/"*.hostfree.fq; do
  s=$(basename "$fq" .hostfree.fq)
  sd="$PER_ASV_DIR/mammal/$s"; mkdir -p "$sd"
  fa="$sd/${s}.fa"; dr="$sd/${s}.derep.fa"; un="$sd/${s}.unoise.fa"

  vsearch --fastq_filter "$fq" --fastaout "$fa" --threads $THREADS
  vsearch --derep_fulllength "$fa" --sizeout --relabel "${s}_" --output "$dr" --threads $THREADS
  vsearch --cluster_unoise "$dr" --minsize 5 --unoise_alpha 2 --centroids "$un" --threads $THREADS

  echo -e "$s\t$(grep -c '^>' "$fa" || echo 0)\t$(grep -c '^>' "$dr" || echo 0)\t$(grep -c '^>' "$un" || echo 0)" \
    >> "$REPORTS_DIR/step9_per_sample_asv_mammal.tsv"
done
tsv_total_row "$REPORTS_DIR/step9_per_sample_asv_mammal.tsv"

echo -e "Sample\tFASTAseqs\tDerepSeqs\tUNOISE_ASVs" > "$REPORTS_DIR/step9_per_sample_asv_12S.tsv"
for fq in "$HOSTFREE_ROOT/12S/"*.hostfree.fq; do
  s=$(basename "$fq" .hostfree.fq)
  sd="$PER_ASV_DIR/12S/$s"; mkdir -p "$sd"
  fa="$sd/${s}.fa"; dr="$sd/${s}.derep.fa"; un="$sd/${s}.unoise.fa"

  vsearch --fastq_filter "$fq" --fastaout "$fa" --threads $THREADS
  vsearch --derep_fulllength "$fa" --sizeout --relabel "${s}_" --output "$dr" --threads $THREADS
  vsearch --cluster_unoise "$dr" --minsize 5 --unoise_alpha 2 --centroids "$un" --threads $THREADS

  echo -e "$s\t$(grep -c '^>' "$fa" || echo 0)\t$(grep -c '^>' "$dr" || echo 0)\t$(grep -c '^>' "$un" || echo 0)" \
    >> "$REPORTS_DIR/step9_per_sample_asv_12S.tsv"
done
tsv_total_row "$REPORTS_DIR/step9_per_sample_asv_12S.tsv"

# =========================
# STEP 10 — MERGE LAST: union per-sample UNOISE → derep → global chimera
# =========================
# mammal
cat "$PER_ASV_DIR"/mammal/*/*.unoise.fa > "$GLOBAL_DIR/mammal/all_unoise.fa" || : 
vsearch --derep_fulllength "$GLOBAL_DIR/mammal/all_unoise.fa" \
        --sizeout --relabel "GMiM_" --output "$GLOBAL_DIR/mammal/all_unoise_uniq.fa" --threads $THREADS || :
vsearch --uchime3_denovo "$GLOBAL_DIR/mammal/all_unoise_uniq.fa" \
        --nonchimeras "$CHIMERA_DIR/ASVs_MiMammal_nochimeras.fasta" --threads $THREADS || :

# 12S
cat "$PER_ASV_DIR"/12S/*/*.unoise.fa > "$GLOBAL_DIR/12S/all_unoise.fa" || :
vsearch --derep_fulllength "$GLOBAL_DIR/12S/all_unoise.fa" \
        --sizeout --relabel "G12S_" --output "$GLOBAL_DIR/12S/all_unoise_uniq.fa" --threads $THREADS || :
vsearch --uchime3_denovo "$GLOBAL_DIR/12S/all_unoise_uniq.fa" \
        --nonchimeras "$CHIMERA_DIR/ASVs_12SV5_nochimeras.fasta" --threads $THREADS || :

echo -e "Locus\tPerSampleUNOISE_Total\tGlobal_Uniq\tGlobal_NoChimera" > "$REPORTS_DIR/step10_global_counts.tsv"
echo -e "MiMammal\t$(grep -ch '^>' "$PER_ASV_DIR"/mammal/*/*.unoise.fa 2>/dev/null | awk '{s+=$1}END{print s+0}')\t$(grep -c '^>' "$GLOBAL_DIR/mammal/all_unoise_uniq.fa" 2>/dev/null || echo 0)\t$(grep -c '^>' "$CHIMERA_DIR/ASVs_MiMammal_nochimeras.fasta" 2>/dev/null || echo 0)" >> "$REPORTS_DIR/step10_global_counts.tsv"
echo -e "12SV5\t$(grep -ch '^>' "$PER_ASV_DIR"/12S/*/*.unoise.fa   2>/dev/null | awk '{s+=$1}END{print s+0}')\t$(grep -c '^>' "$GLOBAL_DIR/12S/all_unoise_uniq.fa" 2>/dev/null || echo 0)\t$(grep -c '^>' "$CHIMERA_DIR/ASVs_12SV5_nochimeras.fasta" 2>/dev/null || echo 0)" >> "$REPORTS_DIR/step10_global_counts.tsv"

# =========================
# STEP 11 — Taxonomy (SINTAX on MIDORI UNIQ)
# =========================
if [[ ! -s "$MIDORI_UDB" ]]; then
  need "$MIDORI_FASTA"
  vsearch --makeudb_usearch "$MIDORI_FASTA" --output "$MIDORI_UDB"
fi

vsearch -sintax "$CHIMERA_DIR/ASVs_MiMammal_nochimeras.fasta" \
        --db "$MIDORI_UDB" --strand both --dbmask none \
        --sintax_cutoff $SINTAX_CUTOFF \
        --tabbedout "$TAX_DIR/ASVs_MiMammal_sintax.tsv" \
        --threads $THREADS || :

# 12S
vsearch -sintax "$CHIMERA_DIR/ASVs_12SV5_nochimeras.fasta" \
        --db "$MIDORI_UDB" --strand both --dbmask none \
        --sintax_cutoff $SINTAX_CUTOFF \
        --tabbedout "$TAX_DIR/ASVs_12SV5_sintax.tsv" \
        --threads $THREADS || :
# =========================
# STEP 12 — Map & build final OTU tables
# =========================

# 12A) Robust ID maps & OTU FASTAs: create valid EMPTY placeholders if needed
make_otu_map () {
  local chim="$1" map="$2" fa="$3"
  mkdir -p "$(dirname "$map")" "$(dirname "$fa")"
  if [[ -s "$chim" ]] && grep -q '>' "$chim"; then
    : > "$map"
    awk 'BEGIN{OFS="\t"}
         /^>/{
             h=substr($0,2);
             split(h,a,/[[:space:];]/); old=a[1];
             n++; print "OTU_" n, old > map; print ">OTU_" n; next
         }
         {print}
    ' map="$map" "$chim" > "$fa"
  else
    : > "$map"; : > "$fa"
  fi
}

M_OTU_MAP="$GLOBAL_DIR/mammal/otu_id_map_mammal.tsv"
M_OTU_FA="$RESULTS_DIR/mammal_otu.fa"
S_OTU_MAP="$GLOBAL_DIR/12S/otu_id_map_12S.tsv"
S_OTU_FA="$RESULTS_DIR/12S_otu.fa"

make_otu_map "$CHIMERA_DIR/ASVs_MiMammal_nochimeras.fasta" "$M_OTU_MAP" "$M_OTU_FA"
make_otu_map "$CHIMERA_DIR/ASVs_12SV5_nochimeras.fasta"   "$S_OTU_MAP" "$S_OTU_FA"

# 12B) Remap SINTAX to new OTU IDs (normalize IDs)
awk 'BEGIN{FS=OFS="\t"} FNR==NR{m[$2]=$1; next}
     {k=$1; sub(/[ ;].*$/,"",k); $1=(k in m? m[k]:$1); print}' \
    "$M_OTU_MAP" "$TAX_DIR/ASVs_MiMammal_sintax.tsv" > "$RESULTS_DIR/mammal_taxa.tsv" || : 

awk 'BEGIN{FS=OFS="\t"} FNR==NR{m[$2]=$1; next}
     {k=$1; sub(/[ ;].*$/,"",k); $1=(k in m? m[k]:$1); print}' \
    "$S_OTU_MAP" "$TAX_DIR/ASVs_12SV5_sintax.tsv"  > "$RESULTS_DIR/12S_taxa.tsv" || :

# =========================
# STEP 12B.1 — Rebuild lineage from full taxonomy (DROP probs & numeric IDs; ':'→'_')
# =========================
rebuild_lineage_from_tax () {
    in_tsv="$1"
    cutoff="${2:-${SINTAX_LIN_CUTOFF:-$SINTAX_CUTOFF}}"
    mode="${3:-${SINTAX_LIN_MODE:-truncate}}"        # optional arg: truncate|filter
    tmp="${in_tsv}.tmp"

    awk -v cutoff="$cutoff" -v mode="$mode" '
        BEGIN{FS=OFS="\t"}
        {
        id=$1
        tax=$2                       # full SINTAX taxonomy (unchanged)
        plus=(NF>=3 ? $3 : "+")      # keep literal "+" if present; insert if missing

      # Build confidence-aware lineage:
      lin=""
      n=split(tax, tok, / *, */)   # split on commas, allow spaces
      for(i=1;i<=n;i++){
        r=""; name=""; conf=-1
        # token looks like: k:Eukaryota_2759(1.00)  or  s:Steno(0.28)
        if (match(tok[i], /^([A-Za-z]):([^,(]+)(\([0-9.]+\))?/, m)) {
          r=m[1]
          name=m[2]
          sub(/_[0-9]+$/, "", name)          # drop numeric IDs like _2759
          gsub(/[[:space:]]+/, "", name)     # be safe about stray spaces
          if (match(tok[i], /\(([0-9.]+)\)/, c)) { conf=c[1]+0 } else { conf=-1 }

          if (conf>=cutoff || conf<0) {
            lin = (lin=="" ? r"_"name : lin","r"_"name)
          } else if (mode=="truncate") {
            break
          } else {
            # mode == "filter": just skip this rank and continue
          }
        }
      }

      print id, tax, plus, lin
    }
  ' "$in_tsv" > "$tmp" && mv -f "$tmp" "$in_tsv"
}


rebuild_lineage_from_tax "$RESULTS_DIR/mammal_taxa.tsv"
rebuild_lineage_from_tax "$RESULTS_DIR/12S_taxa.tsv"

# 12C) Per-sample read mapping (0.99)
for fq in "$HOSTFREE_ROOT/mammal/"*.hostfree.fq; do
  s=$(basename "$fq" .hostfree.fq)
  vsearch --usearch_global "$fq" -db "$CHIMERA_DIR/ASVs_MiMammal_nochimeras.fasta" \
          -id 0.99 -otutabout "$MAP_DIR/mammal/${s}.otutab.tsv" --threads $THREADS || : 
done
for fq in "$HOSTFREE_ROOT/12S/"*.hostfree.fq; do
  s=$(basename "$fq" .hostfree.fq)
  vsearch --usearch_global "$fq" -db "$CHIMERA_DIR/ASVs_12SV5_nochimeras.fasta" \
          -id 0.99 -otutabout "$MAP_DIR/12S/${s}.otutab.tsv" --threads $THREADS || :
done

# 12D) Build OTU tables (resilient; never drops samples due to missing map)
python3 - <<'PY'
import os, csv, glob
from collections import OrderedDict

RESULTS_DIR = os.environ["RESULTS_DIR"]
MAP_DIR     = os.environ["MAP_DIR"]
GLOBAL_DIR  = os.environ["GLOBAL_DIR"]

def norm_id(x):
    if not x: return x
    x = x.split()[0]
    x = x.split(";")[0]
    return x

def write_empty(out_path, sample_tabs, header_label="#OTU ID"):
    samples=[]
    for tsv in sorted(sample_tabs):
        s=os.path.basename(tsv).replace(".otutab.tsv","")
        samples.append(s)
    samples=sorted(list(OrderedDict.fromkeys(samples)))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path,"w",newline="") as f:
        w=csv.writer(f, delimiter="\t"); w.writerow([header_label]+samples)

def build_table(sample_tabs, idmap_path, out_path):
    if not os.path.exists(idmap_path) or os.path.getsize(idmap_path)==0:
        write_empty(out_path, sample_tabs); return

    old2new={}
    with open(idmap_path, newline="") as f:
        for new, old in csv.reader(f, delimiter="\t"):
            old2new[norm_id(old)] = new
    master=[row[0] for row in csv.reader(open(idmap_path), delimiter="\t")]
    table={k:{} for k in master}
    samples=[]
    for tsv in sorted(sample_tabs):
        s=os.path.basename(tsv).replace(".otutab.tsv","")
        samples.append(s)
        if not os.path.getsize(tsv): continue
        with open(tsv, newline="") as f:
            r=csv.reader(f, delimiter="\t"); next(r,None)
            for row in r:
                if not row or len(row)<2: continue
                try: cnt=int(float(row[1]))
                except: cnt=0
                if cnt<=0: continue
                new=old2new.get(norm_id(row[0]))
                if new: table[new][s]=table[new].get(s,0)+cnt
    samples=sorted(list(OrderedDict.fromkeys(samples)))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path,"w",newline="") as f:
        w=csv.writer(f, delimiter="\t")
        w.writerow(["#OTU ID"]+samples)
        for mid in master:
            w.writerow([mid]+[table[mid].get(s,0) for s in samples])

build_table(sorted(glob.glob(os.path.join(MAP_DIR,"mammal","*.otutab.tsv"))),
            os.path.join(GLOBAL_DIR,"mammal","otu_id_map_mammal.tsv"),
            os.path.join(RESULTS_DIR,"mammal_otu_table.tsv"))

build_table(sorted(glob.glob(os.path.join(MAP_DIR,"12S","*.otutab.tsv"))),
            os.path.join(GLOBAL_DIR,"12S","otu_id_map_12S.tsv"),
            os.path.join(RESULTS_DIR,"12S_otu_table.tsv"))
PY


# -------------------------
# NEW: GLOBAL sequence-based collapse (100% identity)
#       Remove Homo + filter per-sample rel. abundance >= 0.1% (0.001)
#       Uses global OTU FASTA + all per-sample *.otutab.tsv
#       Produces CONSOTU_* tables/taxonomy/fasta/map
# -------------------------
# =========================
# STEP 13 — Read-tracking summary (no host columns)
# =========================
python3 - <<'PY'
import csv, os
RESULTS=os.environ["RESULTS_DIR"]
REPORTS=os.environ["REPORTS_DIR"]

def load_simple(p):
    d={}
    if not os.path.exists(p): return d
    with open(p, newline="") as f:
        r=csv.reader(f, delimiter="\t"); next(r,None)
        for row in r:
            if not row or row[0]=="Total": continue
            try: d[row[0]]=int(float(row[1]))
            except: d[row[0]]=0
    return d
def col_sums(p):
    sums={}
    if not os.path.exists(p): return sums
    with open(p, newline="") as f:
        r=csv.reader(f, delimiter="\t"); header=next(r,None)
        if not header: return sums
        samples=header[1:]; sums={s:0 for s in samples}
        for row in r:
            for i,s in enumerate(samples,1):
                try: sums[s]+=int(float(row[i]))
                except: sums[s]+=0
    return sums

m_merged = load_simple(os.path.join(REPORTS,"step6_merged_counts_mammal.tsv"))
s_merged = load_simple(os.path.join(REPORTS,"step6_merged_counts_12S.tsv"))
m_map    = col_sums   (os.path.join(RESULTS,"mammal_otu_table.tsv"))
s_map    = col_sums   (os.path.join(RESULTS,"12S_otu_table.tsv"))

samples = sorted(set(m_merged)|set(s_merged)|set(m_map)|set(s_map))
out = os.path.join(RESULTS,"read_tracking_summary.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out,"w",newline="") as f:
    w=csv.writer(f, delimiter="\t")
    w.writerow(["Sample","Merged_mammal","Mapped_mammal",
                "Merged_12S","Mapped_12S",
                "Total_Merged","Total_Mapped"])
    for s in samples:
        mm=m_merged.get(s,0); mx=m_map.get(s,0)
        sm=s_merged.get(s,0); sx=s_map.get(s,0)
        w.writerow([s,mm,mx,sm,sx,mm+sm,mx+sx])
PY

# =========================
# STEP 12V/13.5 — Integrity checks (compare OTU table sums vs tracking)
# =========================
audit_locus () {
  local locus="$1"
  if [[ "$locus" == "12S" ]]; then
    local CHIM="$CHIMERA_DIR/ASVs_12SV5_nochimeras.fasta"
    local MAP="$GLOBAL_DIR/12S/otu_id_map_12S.tsv"
    local OTUTAB="$RESULTS_DIR/12S_otu_table.tsv"
    local TAXA="$RESULTS_DIR/12S_taxa.tsv"
    local UNOISE="$PER_ASV_DIR/12S/*/*.unoise.fa"
    local ALLUNIQ="$GLOBAL_DIR/12S/all_unoise_uniq.fa"
    local MAPPED_COL="Mapped_12S"
  else
    local CHIM="$CHIMERA_DIR/ASVs_MiMammal_nochimeras.fasta"
    local MAP="$GLOBAL_DIR/mammal/otu_id_map_mammal.tsv"
    local OTUTAB="$RESULTS_DIR/mammal_otu_table.tsv"
    local TAXA="$RESULTS_DIR/mammal_taxa.tsv"
    local UNOISE="$PER_ASV_DIR/mammal/*/*.unoise.fa"
    local ALLUNIQ="$GLOBAL_DIR/mammal/all_unoise_uniq.fa"
    local MAPPED_COL="Mapped_mammal"
  fi

  local UNOISE_TOTAL=$(grep -ch '^>' $UNOISE 2>/dev/null | awk '{s+=$1}END{print s+0}')
  local GLOBAL_UNIQ=$(grep -c '^>' "$ALLUNIQ" 2>/dev/null || echo 0)
  local NOCHIM=$(grep -c '^>' "$CHIM" 2>/dev/null || echo 0)
  local MAP_LINES=$(wc -l < "$MAP" 2>/dev/null || echo 0)
  local OTU_ROWS=$(( $(wc -l < "$OTUTAB" 2>/dev/null || echo 1) - 1 ))
  local TAXA_ROWS=$(wc -l < "$TAXA" 2>/dev/null || echo 0)

  {
    echo -e "Metric\tCount"
    echo -e "PerSample_UNOISE_total\t$UNOISE_TOTAL"
    echo -e "Global_unique_after_derep\t$GLOBAL_UNIQ"
    echo -e "Nonchimeric_ASVs\t$NOCHIM"
    echo -e "OTU_map_rows\t$MAP_LINES"
    echo -e "OTU_table_rows\t$OTU_ROWS"
    echo -e "Taxa_rows\t$TAXA_ROWS"
    [[ "$MAP_LINES" -eq "$NOCHIM" ]] || echo -e "WARNING_map_vs_nochim\t${MAP_LINES}!="${NOCHIM}
    [[ "$OTU_ROWS" -eq "$MAP_LINES" ]] || echo -e "WARNING_table_vs_map\t${OTU_ROWS}!="${MAP_LINES}
    [[ "$TAXA_ROWS" -eq "$MAP_LINES" ]] || echo -e "WARNING_taxa_vs_map\t${TAXA_ROWS}!="${MAP_LINES}
  } > "$REPORTS_DIR/step12_integrity_${locus}.tsv"

  awk -F'\t' 'NR==1{for(i=2;i<=NF;i++) h[i]=$i; next}
              {for(i=2;i<=NF;i++) s[h[i]]+=$i}
              END{for(k in s) printf "%s\t%d\n", k, s[k]}' "$OTUTAB" \
    | sort -k1,1 > "$REPORTS_DIR/otutab_sums_${locus}.tsv" || :

  awk -F'\t' -v col="$MAPPED_COL" '
    NR==1{for(i=1;i<=NF;i++) if($i==col) C=i; next}
    C && NR>1{printf "%s\t%d\n",$1, ($C==""?0:$C)}
  ' "$RESULTS_DIR/read_tracking_summary.tsv" | sort -k1,1 > "$REPORTS_DIR/tracked_${locus}.tsv" || :

  join -t $'\t' -a1 -a2 -e 0 -o 0,1.2,2.2 \
    "$REPORTS_DIR/otutab_sums_${locus}.tsv" "$REPORTS_DIR/tracked_${locus}.tsv" \
    | awk -F'\t' '$2!=$3{print $0}' > "$REPORTS_DIR/step12_integrity_${locus}_mismatch.tsv" || true
}

audit_locus 12S
audit_locus mammal

# =========================
# apply the filter 
# =========================
python3 $SCRIPTS_DIR/collapse_by_sequence_global_f.py \
  --input-dir $RESULTS_DIR \
  --out-dir $RESULTS_DIR \
  --rel-threshold 0.01 \
  --control-list "NCDNAmerged_mammal,NCPCRmerged_mammal,PCmerged_mammal"

# after your collapse step
python3 "$SCRIPTS_DIR/report.py" --input-dir "$RESULTS_DIR"


# =========================
# DONE
# =========================
echo "🎉 DONE. Final outputs in $RESULTS_DIR :"
ls -lh "$RESULTS_DIR"/{mammal_otu.fa,12S_otu.fa,mammal_otu_table.tsv,12S_otu_table.tsv,mammal_taxa.tsv,12S_taxa.tsv,read_tracking_summary.tsv,mammal_sample_metrics.tsv,12S_sample_metrics.tsv} 2>/dev/null || true
echo "ℹ️ See integrity & QC reports in $REPORTS_DIR:"
ls -lh "$REPORTS_DIR"/step12_integrity_* "$REPORTS_DIR/step14_qc_otu_inflation.tsv" 2>/dev/null || true
