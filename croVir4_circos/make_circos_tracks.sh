#!/usr/bin/env bash
set -euo pipefail

FASTA="croVir4_genome_2025-10-05_v1.fasta"
GTF="croVir4_annotation_2026-03-02_v1.gtf"
WINDOW=500000
OUTDIR="circos_tracks"

mkdir -p "$OUTDIR"

echo "Generating genome windows (${WINDOW} bp)"
samtools faidx "$FASTA"
cut -f1,2 "${FASTA}.fai" > genome.sizes
bedtools makewindows -g genome.sizes -w "$WINDOW" > "${OUTDIR}/windows.bed"

echo "Calculating GC content"
bedtools nuc -fi "$FASTA" -bed "${OUTDIR}/windows.bed" \
  | awk 'NR>1 {printf "%s\t%s\t%s\t%.4f\n", $1, $2, $3, $5}' \
  > "${OUTDIR}/gc_content.txt"

echo "Extracting gene features from GTF (using transcripts as proxy)"
awk '$3 == "transcript" {OFS="\t"; print $1, $4-1, $5}' "$GTF" \
  | sort -k1,1 -k2,2n \
  > "${OUTDIR}/genes.bed"

echo "==> Calculating gene density"
bedtools coverage -a "${OUTDIR}/windows.bed" -b "${OUTDIR}/genes.bed" \
  | awk '{print $1, $2, $3, $4}' OFS="\t" \
  > "${OUTDIR}/gene_density.txt"

echo "==> Calculating repeat content"
REPEAT_GFF="croVir4_genome_2026-04-29_v1_RepeatFamily.gff3"
grep -v '^#' "$REPEAT_GFF" \
  | awk '{OFS="\t"; print $1, $4-1, $5}' \
  | sort -k1,1 -k2,2n \
  > "${OUTDIR}/repeats.bed"
bedtools coverage -a "${OUTDIR}/windows.bed" -b "${OUTDIR}/repeats.bed" \
  | awk '{printf "%s\t%s\t%s\t%.4f\n", $1, $2, $3, $7}' \
  > "${OUTDIR}/repeat_content.txt"

rm -f genome.sizes
