#!/bin/bash
# CPU vs GPU alignment benchmark on a few samples.
#
#   1. fastp once per sample; every arm aligns the same trimmed reads, and
#      trimming is not part of the timed comparison.
#   2. CPU arms (STAR binary + matching index), CPU_PAR samples at a time x
#      CPU_THREADS threads.
#   3. GPU arms: Parabricks rna_fq2bam on n GPUs (GPU_COUNTS), one sample at a time.
#   4. bin/benchmark_compare.py: times, mapping rates, junction and gene-count agreement.
# All arms use MISMATCH="--outFilterMismatchNoverLmax 0.04" (see bin/star_align.sh),
# the closest setting Parabricks offers to the production filter.
# GPU steps start no earlier than GPU_START (HH:MM, today) and only on GPUs with no
# running processes. Finished samples are skipped, so the script can be rerun.
#
# Usage: benchmark_align.sh <samplesheet.csv> <raw_dir> <index_2.7.2a> <genome.fa>
#                           <outdir> <sample,sample,...> <arm=STAR:index> ...
#   e.g. cpu_2.7.11b=/path/STAR:/path/index_2.7.11b  (index_2.7.2a is the GPU index)
# Env: PARABRICKS_SIF, ADAPTER, CPU_PAR (1), CPU_THREADS (48), GPU_START (19:00),
#      GPU_COUNTS ("1")
set -euo pipefail
SHEET=$1 RAW=$2 IDX272=$3 FASTA=$4 OUT=$5 BENCH=$6; shift 6
CPU_ARMS=("$@")
: "${PARABRICKS_SIF:?}" "${ADAPTER:?}"
CPU_PAR=${CPU_PAR:-1} CPU_THREADS=${CPU_THREADS:-48} GPU_START=${GPU_START:-19:00}
GPU_COUNTS=${GPU_COUNTS-1}   # empty: no GPU arms
BIN=$(cd "$(dirname "$0")" && pwd)
export MISMATCH="--outFilterMismatchNoverLmax 0.04" PARABRICKS_SIF
T=$OUT/trim TIMING=$OUT/timing.tsv
mkdir -p "$T" "$OUT/fastp"
[ -s "$TIMING" ] || printf 'arm\tsample\twall_s\tmax_rss_kb\tstart\tend\n' > "$TIMING"
log() { echo "[$(date '+%F %T')] $*"; }

secs() {  # /usr/bin/time -v "Elapsed (wall clock)" -> seconds
    awk -F': ' '/Elapsed/ { n = split($2, t, ":"); s = 0; for (i = 1; i <= n; i++) s = s * 60 + t[i]; print s }' "$1"
}
rss() { awk -F': ' '/Maximum resident/ { print $2 }' "$1"; }
record() {  # arm sample start end time_file
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$(secs "$5")" "$(rss "$5")" "$3" "$4" >> "$TIMING"
}
done_bam() { [ -s "$1.Aligned.sortedByCoord.out.bam.bai" ] && [ -s "$1.Log.final.out" ]; }
wait_gpus() {  # n -> prints n idle GPU indices, once GPU_START has passed
    local n=$1 idle busy
    while :; do
        if [ "$(date +%H%M)" -ge "${GPU_START/:/}" ]; then
            busy=$(nvidia-smi --query-compute-apps=gpu_uuid --format=csv,noheader | sort -u)
            idle=$(nvidia-smi --query-gpu=index,uuid --format=csv,noheader | sed 's/, /,/' \
                   | while IFS=, read -r i u; do grep -qx "$u" <<< "$busy" || echo "$i"; done \
                   | head -n "$n" | paste -sd,)
            [ "$(tr ',' '\n' <<< "$idle" | grep -c .)" -eq "$n" ] && { echo "$idle"; return; }
        fi
        sleep 60
    done
}
IFS=, read -r -a SAMPLES <<< "$BENCH"

# 1. Trim once.
for s in "${SAMPLES[@]}"; do
    [ -s "$T/${s}_R2.fq.gz" ] && [ -s "$OUT/fastp/$s.fastp.json" ] && continue
    read -r r1 r2 < <(python3 -I -c 'import csv,sys
for r in csv.DictReader(open(sys.argv[1])):
    if r["sample"] == sys.argv[2]: print(r["fastq_1"], r["fastq_2"])' "$SHEET" "$s")
    log "fastp $s"
    /usr/bin/time -v -o "$OUT/fastp/$s.time.txt" fastp -i "$RAW/$r1" -I "$RAW/$r2" \
        -o "$T/${s}_R1.fq.gz" -O "$T/${s}_R2.fq.gz" \
        --adapter_sequence "$ADAPTER" --adapter_sequence_r2 "$ADAPTER" \
        --disable_quality_filtering --length_required 20 --trim_poly_g -w 8 \
        -j "$OUT/fastp/$s.fastp.json" -h "$OUT/fastp/$s.fastp.html" 2> "$OUT/fastp/$s.fastp.log"
done

# 2. CPU arms.
cpu_one() {  # arm star index sample; runs under xargs bash -c, so no set -e: check exits
    local d=$OUT/$1/$4 st en
    done_bam "$d/$4" && return 0
    st=$(date +%s)
    STAR=$2 "$BIN/star_align.sh" "$3" "$4" "$T/${4}_R1.fq.gz" "$T/${4}_R2.fq.gz" "$d" "$CPU_THREADS" \
        > "$OUT/$1.$4.log" 2>&1 || { echo "FAILED $1 $4 (see $OUT/$1.$4.log)" >&2; return 1; }
    en=$(date +%s); record "$1" "$4" "$st" "$en" "$d/$4.time.txt"
}
export -f cpu_one record secs rss done_bam
export OUT T BIN TIMING CPU_THREADS
for a in "${CPU_ARMS[@]}"; do
    arm=${a%%=*} spec=${a#*=}; star=${spec%%:*} idx=${spec#*:}
    log "$arm: ${#SAMPLES[@]} samples, $CPU_PAR x $CPU_THREADS threads"
    printf '%s\n' "${SAMPLES[@]}" | xargs -P "$CPU_PAR" -I{} bash -c 'cpu_one "$@"' _ "$arm" "$star" "$idx" {}
done

# 3. GPU arms.
for n in $GPU_COUNTS; do
    arm=gpu_$n
    for s in "${SAMPLES[@]}"; do
        d=$OUT/$arm/$s; done_bam "$d/$s" && continue
        g=$(wait_gpus "$n"); log "$arm $s on GPU $g"
        st=$(date +%s)
        CUDA_VISIBLE_DEVICES=$g "$BIN/parabricks_align.sh" "$IDX272" "$FASTA" "$s" \
            "$T/${s}_R1.fq.gz" "$T/${s}_R2.fq.gz" "$d" "$n" > "$OUT/$arm.$s.log" 2>&1
        en=$(date +%s); record "$arm" "$s" "$st" "$en" "$d/$s.time.txt"
    done
done

# 4. Compare.
python3 -I "$BIN/benchmark_compare.py" "$OUT" > "$OUT/comparison.txt"
rm -f "$T"/*.fq.gz
log "done"; cat "$OUT/comparison.txt"
