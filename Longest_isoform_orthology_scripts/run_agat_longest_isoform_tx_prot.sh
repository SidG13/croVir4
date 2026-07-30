#!/bin/bash
# Usage: ./run_agat_longest_isoform.sh <input.gtf> <genome.fasta>
# Example: ./run_agat_longest_isoform.sh croVir4_annotation_2026-06-12_v1.gtf genome.fasta

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.gtf> <genome.fasta>"
    exit 1
fi

INPUT_GTF=$1
GENOME_FASTA=$2

# Strip .gtf extension to get the base name
BASE=$(basename "$INPUT_GTF" .gtf)

OUT_GTF="${BASE}_KEEPLONGESTISOFORM-1.gtf"
OUT_PROT="${BASE}_KEEPLONGESTISOFORM-PROTEINS.fasta"
OUT_TRANS="${BASE}_KEEPLONGESTISOFORM-TRANSCIPTS.fasta"

echo "Step 1: Keeping longest isoform..."
agat_sp_keep_longest_isoform.pl \
    --gff "$INPUT_GTF" \
    -o "$OUT_GTF"

echo "Step 2: Extracting protein sequences..."
agat_sp_extract_sequences.pl \
    -g "$OUT_GTF" \
    -f "$GENOME_FASTA" \
    -t cds \
    -p \
    -o "$OUT_PROT"

echo "Step 3: Extracting transcript sequences..."
agat_sp_extract_sequences.pl \
    -g "$OUT_GTF" \
    -f "$GENOME_FASTA" \
    -t exon \
    --merge \
    -o "$OUT_TRANS"

echo "Done. Outputs:"
echo "  $OUT_GTF"
echo "  $OUT_PROT"
echo "  $OUT_TRANS"
