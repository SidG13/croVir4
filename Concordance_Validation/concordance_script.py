import csv
from collections import defaultdict

# --- 1. Parse PAF, keep best (highest-coverage) alignment per contig ---
final_pos = {}
contig_len = {}
with open("draft_to_final.paf") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        contig, qlen, qstart, qend = p[0], int(p[1]), int(p[2]), int(p[3])
        final_chrom = p[5]
        mapq = int(p[11])
        cov = (qend - qstart) / qlen
        if mapq < 20 or cov < 0.5:
            continue
        if contig not in final_pos or cov > final_pos[contig][1]:
            final_pos[contig] = (final_chrom, cov)
        contig_len[contig] = qlen

# --- 2. YaHS scaffold -> ordered contig list ---
yahs_scaffolds = defaultdict(list)
with open("yahs_order.tsv") as f:
    for row in csv.reader(f, delimiter="\t"):
        scaffold, contig, strand = row
        yahs_scaffolds[scaffold].append(contig)

# --- 3. Per-scaffold bp-weighted majority vote ---
results = []
for scaffold, contigs in yahs_scaffolds.items():
    bp_by_chrom = defaultdict(int)
    n = 0
    for c in contigs:
        if c not in final_pos or c not in contig_len:
            continue
        chrom = final_pos[c][0]
        bp_by_chrom[chrom] += contig_len[c]
        n += 1
    if n == 0:
        continue
    majority_chrom = max(bp_by_chrom, key=bp_by_chrom.get)
    total_bp = sum(bp_by_chrom.values())
    pct_concordant = bp_by_chrom[majority_chrom] / total_bp   # 100% = fully concordant
    results.append((scaffold, majority_chrom, pct_concordant, n, len(bp_by_chrom)))

# --- 4. Report, best-to-worst (100% concordant first) ---
results.sort(key=lambda r: -r[2])

print(f"{'scaffold':15s} {'majority_chrom':16s} {'%concordant':>12s} {'n_contigs':>10s} {'n_chroms':>9s}")
for scaffold, majority_chrom, pct, n, n_chroms in results:
    flag = "  <-- LOW n" if n < 10 else ""
    print(f"{scaffold:15s} {majority_chrom:16s} {pct:12.1%} {n:10d} {n_chroms:9d}{flag}")

concordant = [r for r in results if r[2] >= 0.95]
print(f"\n{len(concordant)}/{len(results)} multi-contig scaffolds >=95% concordant")

with open("scaffold_concordance_final.tsv", "w") as out:
    out.write("scaffold\tmajority_chrom\tpct_concordant\tn_contigs\tn_chroms\n")
    for r in results:
        out.write(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]}\t{r[4]}\n")