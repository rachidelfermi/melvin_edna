#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --account=omics_core-rqyo8fdrbpw-LOW-CPU
#SBATCH --qos=low-cpu
#SBATCH --job-name=melvin_demux_merge

# === Load environment ===
eval "$(/srv/data/rachid.elfermi/miniconda3/bin/conda shell.bash hook)"
conda activate bioinfo

cutadapt

## === Setup working directory ===
WORKDIR="/home/rachid.elfermi/lustre/omics_core-rqyo8fdrbpw/users/rachid.elfermi/dada2_nfcore/melvin"
cd "$WORKDIR" || { echo "❌ Failed to enter $WORKDIR"; exit 1; }

# === Create folders ===
mkdir -p 02_qualityfilter_beht 03_mergepairs_beht 04_demultiplex_level1_beht 05_demultiplex_level2_beht
mkdir -p 05_demultiplex_level2_beht/{MiMammal1,MiMammal2,12SV51,12SV52}
mkdir -p 06_merged_beht/mammal 06_merged_beht/12S 07_reports_beht

# === Helpers ===
count_reads() {
    # Count FASTQ records by counting '+' lines (robust to long headers)
    local f="$1"
    if [[ -f "$f" ]]; then
        if [[ "$f" == *.gz ]]; then
            zcat -f -- "$f" 2>/dev/null | grep -c "^+$" || echo 0
        else
            grep -c "^+$" -- "$f" 2>/dev/null || echo 0
        fi
    else
        echo 0
    fi
}

# Input filenames (adjust here if names change)
R1_IN="Melvin_FKDN240276974-1A_HLY2WDSXC_L2_1.fq"
R2_IN="Melvin_FKDN240276974-1A_HLY2WDSXC_L2_2.fq"

R1_FILT="02_qualityfilter_beht/Melvin_1.fq"
R2_FILT="02_qualityfilter_beht/Melvin_2.fq"
MERGED="03_mergepairs_beht/Melvin_merged.fq"

# === (Preview) Second report — original counts before any processing ===
# We’ll fill this again after each stage to reflect actual outputs
echo -e "original_R1\toriginal_R2\tfiltered_R1\tfiltered_R2\tmerged_pairs" > 07_reports_beht/run_overview.tsv
orig_r1=$(count_reads "$R1_IN")
orig_r2=$(count_reads "$R2_IN")
echo -e "${orig_r1}\t${orig_r2}\t0\t0\t0" >> 07_reports_beht/run_overview.tsv

# === Step 1: Quality filtering ===
vsearch --fastx_filter "$R1_IN" \
        --reverse "$R2_IN" \
        --fastq_maxee 1 \
        --fastqout "$R1_FILT" \
        --fastqout_rev "$R2_FILT" \
        --threads 10

# Update second report after filtering
filt_r1=$(count_reads "$R1_FILT")
filt_r2=$(count_reads "$R2_FILT")
# overwrite the second line with updated values keeping original_* from first line
# (simpler: rewrite the whole file)
echo -e "original_R1\toriginal_R2\tfiltered_R1\tfiltered_R2\tmerged_pairs" > 07_reports_beht/run_overview.tsv
echo -e "${orig_r1}\t${orig_r2}\t${filt_r1}\t${filt_r2}\t0" >> 07_reports_beht/run_overview.tsv

# === Step 2: Merge pairs ===
vsearch --fastq_mergepairs "$R1_FILT" \
        --reverse "$R2_FILT" \
        --fastqout "$MERGED" \
        --fastq_allowmergestagger \
        --threads 10

# Update second report after merging
merged_pairs=$(count_reads "$MERGED")
echo -e "original_R1\toriginal_R2\tfiltered_R1\tfiltered_R2\tmerged_pairs" > 07_reports_beht/run_overview.tsv
echo -e "${orig_r1}\t${orig_r2}\t${filt_r1}\t${filt_r2}\t${merged_pairs}" >> 07_reports_beht/run_overview.tsv

# === Step 3: Demultiplexing Level 1 ===
cutadapt -g file:./Beht_tags_level1_bothways.fasta \
         -o "04_demultiplex_level1_beht/{name}.fastq" \
         "$MERGED" \
         -j 16 -e 0

# === Step 4: Demultiplexing Level 2 ===
LIB="Melvin"
demul_samples_var="demul_samples_${LIB}"
eval "${demul_samples_var}=\"\$(ls 04_demultiplex_level1_beht | sort | uniq | cut -d. -f1)\""

eval "for s in \$$demul_samples_var; do
    cutadapt -g file:tags_level2_bothways_orientation.fasta \
             -o \"05_demultiplex_level2_beht/\${s}{name}.fastq\" \
             04_demultiplex_level1_beht/\${s}.fastq \
             -j 16 -e 0.125
done"

# === Remove unknowns and files smaller than 1MB ===
echo "🧹 Cleaning demuxed files..."
#rm -f 05_demultiplex_level2_beht/*unknown*.fastq
find 05_demultiplex_level2_beht -type f -name "*.fastq" -size -1M -delete

# === Organize demuxed reads ===
mv 05_demultiplex_level2_beht/*A1.fastq 05_demultiplex_level2_beht/MiMammal1/ 2>/dev/null || true
mv 05_demultiplex_level2_beht/*A2.fastq 05_demultiplex_level2_beht/MiMammal2/ 2>/dev/null || true
mv 05_demultiplex_level2_beht/*B1.fastq 05_demultiplex_level2_beht/12SV51/ 2>/dev/null || true
mv 05_demultiplex_level2_beht/*B2.fastq 05_demultiplex_level2_beht/12SV52/ 2>/dev/null || true

# === Step 5: Reverse complement A2 & B2 ===
echo "🔁 Reverse complementing A2 and B2 FASTQs..."
shopt -s nullglob
for f in 05_demultiplex_level2_beht/MiMammal2/*A2.fastq; do
    base=$(basename "$f" A2.fastq)
    vsearch --fastx_revcomp "$f" \
            --fastqout "05_demultiplex_level2_beht/MiMammal2/${base}A2_RC.fastq"
done

for f in 05_demultiplex_level2_beht/12SV52/*B2.fastq; do
    base=$(basename "$f" B2.fastq)
    vsearch --fastx_revcomp "$f" \
            --fastqout "05_demultiplex_level2_beht/12SV52/${base}B2_RC.fastq"
done
shopt -u nullglob

# === Step 6: Merge MiMammal and 12S reads and track ===
# Now includes demux_level1 column
echo -e "Sample\tdemux_level1\tmammal_f\tmammal_r\tmammal_rc\tmerged_mammal\t12S_f\t12S_r\t12S_rc\tmerged_12S" > 07_reports_beht/read_count_summary.tsv

declare -A counts

# Prefill with demux_level1 counts
for f in 04_demultiplex_level1_beht/*.fastq; do
    base=$(basename "$f" .fastq)
    dmx=$(count_reads "$f")
    counts[$base]="$dmx" # start row with demux level1
done

# Collect mammal reads
shopt -s nullglob
for f in 05_demultiplex_level2_beht/MiMammal1/*A1.fastq; do
    base=$(basename "$f" A1.fastq)
    A1="05_demultiplex_level2_beht/MiMammal1/${base}A1.fastq"
    A2="05_demultiplex_level2_beht/MiMammal2/${base}A2.fastq"
    A2_RC="05_demultiplex_level2_beht/MiMammal2/${base}A2_RC.fastq"
    merged="06_merged_beht/mammal/${base}merged_mammal.fastq"

    reads_a1=$(count_reads "$A1")
    reads_a2=$(count_reads "$A2")
    reads_a2rc=$(count_reads "$A2_RC")

    if [[ -f "$A2_RC" ]]; then
        cat "$A1" "$A2_RC" > "$merged"
        reads_merged_mammal=$(count_reads "$merged")
    else
        reads_merged_mammal=0
    fi

    if [[ -v counts[$base] ]]; then
        counts[$base]="${counts[$base]}	${reads_a1}	${reads_a2}	${reads_a2rc}	${reads_merged_mammal}"
    else
        counts[$base]="0	${reads_a1}	${reads_a2}	${reads_a2rc}	${reads_merged_mammal}"
    fi
done
shopt -u nullglob

# Add 12S reads
shopt -s nullglob
for f in 05_demultiplex_level2_beht/12SV51/*B1.fastq; do
    base=$(basename "$f" B1.fastq)
    B1="05_demultiplex_level2_beht/12SV51/${base}B1.fastq"
    B2="05_demultiplex_level2_beht/12SV52/${base}B2.fastq"
    B2_RC="05_demultiplex_level2_beht/12SV52/${base}B2_RC.fastq"
    merged="06_merged_beht/12S/${base}merged_12S.fastq"

    reads_b1=$(count_reads "$B1")
    reads_b2=$(count_reads "$B2")
    reads_b2rc=$(count_reads "$B2_RC")

    if [[ -f "$B2_RC" ]]; then
        cat "$B1" "$B2_RC" > "$merged"
        reads_merged12s=$(count_reads "$merged")
    else
        reads_merged12s=0
    fi

    if [[ -v counts[$base] ]]; then
        # If mammal section already appended 4 fields, append the 4 for 12S
        counts[$base]="${counts[$base]}	${reads_b1}	${reads_b2}	${reads_b2rc}	${reads_merged12s}"
    else
        # No mammal; ensure we still have demux_level1 field (0 if missing)
        counts[$base]="0	0	0	0	0	${reads_b1}	${reads_b2}	${reads_b2rc}	${reads_merged12s}"
    fi
done
shopt -u nullglob

# Output all samples to summary file
for key in "${!counts[@]}"; do
    echo -e "$key\t${counts[$key]}" >> 07_reports_beht/read_count_summary.tsv
done

# === 📊 Add TOTAL row ===
echo "➕ Adding total read counts..."
awk 'BEGIN{OFS="\t"} NR==1{print; next} $1!="Total"{for(i=2;i<=NF;i++) sum[i]+=$i} END{printf "Total"; for(i=2;i<=NF;i++) printf "\t%d", sum[i]; print ""}' 07_reports_beht/read_count_summary.tsv >> 07_reports_beht/read_count_summary.tsv

echo "✅ All steps complete!"
echo "📊 Per-run overview: 07_reports_beht/run_overview.tsv"
echo "📊 Per-sample tracking: 07_reports_beht/read_count_summary.tsv"
