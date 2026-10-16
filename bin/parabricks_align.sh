#!/bin/bash
# Align one paired-end sample with NVIDIA Parabricks rna_fq2bam (GPU STAR,
# built on STAR 2.7.2a), timed with /usr/bin/time -v. Same settings as
# bin/star_align.sh, option by option, except:
#  - the mismatch cap is relative to the mapped length (Parabricks has no
#    --outFilterMismatchNoverReadLmax); the CPU benchmark arms use the same cap;
#  - duplicates are not marked (--no-markdups), as in the CPU run;
#  - the read group also has PU (platform unit), which Parabricks requires.
# Outputs use the star_align.sh names, so the workflow's --aligned_dir reads them.
#
# Usage: parabricks_align.sh <index_dir> <genome.fa> <sample> <R1.fq.gz> <R2.fq.gz> <out_dir> [n_gpus]
# GPUs: set CUDA_VISIBLE_DEVICES. The container is $PARABRICKS_SIF.
set -euo pipefail
IDX=$1 FASTA=$2 SAMPLE=$3 R1=$4 R2=$5 OUT=$6 NGPU=${7:-1}
SIF=${PARABRICKS_SIF:?set PARABRICKS_SIF to the Parabricks container}
mkdir -p "$OUT" "$OUT/tmp"
apptainer exec "$SIF" pbrun version 2>/dev/null | grep pbrun > "$OUT/star_version.txt"

# Inputs and outputs are bound into the container at the same paths.
binds=$(for p in "$IDX" "$FASTA" "$R1" "$R2" "$OUT"; do readlink -f "$p"; done | cut -d/ -f2 | sort -u | sed 's,.*,/&:/&,' | paste -sd,)

/usr/bin/time -v -o "$OUT/${SAMPLE}.time.txt" \
apptainer exec --nv --bind "$binds" "$SIF" pbrun rna_fq2bam \
  --ref "$(readlink -f "$FASTA")" \
  --genome-lib-dir "$(readlink -f "$IDX")" \
  --in-fq "$(readlink -f "$R1")" "$(readlink -f "$R2")" \
      "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA\tPU:${SAMPLE}" \
  --read-files-command zcat \
  --output-dir "$(readlink -f "$OUT")" --out-prefix "${SAMPLE}." \
  --out-bam "$(readlink -f "$OUT")/${SAMPLE}.Aligned.sortedByCoord.out.bam" \
  --tmp-dir "$(readlink -f "$OUT")/tmp" \
  --num-gpus "$NGPU" --gpusort --gpuwrite \
  --no-markdups \
  --two-pass-mode Basic \
  --out-sam-attributes NH HI AS nM MD \
  --out-sam-strand-field intronMotif \
  --quantMode GeneCounts \
  --max-out-filter-multimap 20 \
  --max-out-filter-mismatch 999 --max-out-filter-mismatch-ratio 0.04 \
  --min-align-sj-overhang 8 --min-align-sjdb-overhang 1 \
  --min-intron-size 20 --max-intron-size 1000000 --max-align-mates-gap 1000000 \
  --logfile "$OUT/${SAMPLE}.parabricks.log"

samtools index -@ 4 "$OUT/${SAMPLE}.Aligned.sortedByCoord.out.bam"
rm -rf "$OUT/tmp" "$OUT/${SAMPLE}._STARgenome" "$OUT/${SAMPLE}._STARpass1" "$OUT/${SAMPLE}._STARtmp"
grep -E 'Elapsed|Maximum resident' "$OUT/${SAMPLE}.time.txt"
