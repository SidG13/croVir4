library(rtracklayer)
library(Biostrings)
library(GenomicRanges)
library(dplyr)

setwd('/Users/sidgopalan/Documents/croVir4_matrix_builder')

gtf <- rtracklayer::import("croVir4_annotation_2026-06-12_v1.gtf")
gtf_df <- as.data.frame(gtf)
venom_gtf <- gtf_df %>% 
  filter(grepl('Venom', gene_name))

genes <- venom_gtf %>% filter(type == "transcript")
exons  <- venom_gtf %>% filter(type == "exon")

exon_summary <- exons %>%
  group_by(gene_name) %>%
  summarise(
    n_exons = n(),
    .groups = "drop"
  )

coords <- genes %>%
  select(gene_name, seqnames, start, end, strand) %>%
  left_join(exon_summary, by = "gene_name") %>%
  mutate(family = sub("^Venom_", "", gene_name) %>% sub("-.*$", "", .))

prot <- readAAStringSet("venom_proteins.fasta")
prot_seq_df <- data.frame(
  gene_name = sub(" .*", "", names(prot)),
  protein_seq = as.character(prot),
  protein_length = width(prot)
) %>%
  mutate(
    starts_with_M = substr(protein_seq, 1, 1) == "M",
    ends_with_stop = substr(protein_seq, nchar(protein_seq), nchar(protein_seq)) == "*",
    has_internal_stop = grepl("\\*", substr(protein_seq, 1, nchar(protein_seq) - 1)),
    cds_bp = (protein_length + 1) * 3
  )

evidence_matrix <- coords %>%
  left_join(prot_seq_df, by = "gene_name")

# SignalP
signalp <- read.delim("output_protein_type.txt", header = FALSE, comment.char = "#",
                      col.names = c("gene_name", "sp_prediction", "sp_prob", "other_prob", "cs_info"))

evidence_matrix <- evidence_matrix %>%
  left_join(signalp %>% select(gene_name, sp_prediction, sp_prob, cs_info), by = "gene_name") %>%
  mutate(has_signal_peptide = sp_prediction == "SP(Sec/SPI)")

# Iso-Seq support (coordinate overlap)
isoseq_gtf <- import("CV1087.Viridis.ISOSEQ.collapsed.croVir4.gtf")
isoseq_tx <- isoseq_gtf[isoseq_gtf$type == "transcript"]

gene_gr <- GRanges(evidence_matrix$seqnames,
                   IRanges(evidence_matrix$start, evidence_matrix$end),
                   strand = evidence_matrix$strand)

hits <- findOverlaps(gene_gr, isoseq_tx)

isoseq_counts <- data.frame(row = queryHits(hits)) %>%
  count(row, name = "n_isoseq_transcripts")

evidence_matrix$n_isoseq_transcripts <- 0L
evidence_matrix$n_isoseq_transcripts[isoseq_counts$row] <- isoseq_counts$n_isoseq_transcripts
evidence_matrix <- evidence_matrix %>%
  mutate(isoseq_supported = n_isoseq_transcripts > 0)
head(evidence_matrix)


# Add information from croVir3 diamond comparison
new_vs_old <- read.delim("croVir4_vs_croVir3_full.tsv", header = FALSE,
                         col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen",
                                       "qstart","qend","sstart","send","evalue","bitscore","qlen","slen"))

old_vs_new <- read.delim("croVir3_vs_croVir4_full.tsv", header = FALSE,
                         col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen",
                                       "qstart","qend","sstart","send","evalue","bitscore","qlen","slen"))

best_new_vs_old <- new_vs_old %>%
  group_by(qseqid) %>% slice_max(bitscore, n = 1, with_ties = FALSE) %>% ungroup()

best_old_vs_new <- old_vs_new %>%
  group_by(qseqid) %>% slice_max(bitscore, n = 1, with_ties = FALSE) %>% ungroup()

rbh_table <- best_new_vs_old %>%
  left_join(best_old_vs_new %>% select(old_gene = qseqid, new_best_hit = sseqid),
            by = c("sseqid" = "old_gene")) %>%
  mutate(
    is_rbh = qseqid == new_best_hit,
    query_coverage = (qend - qstart + 1) / qlen,
    subject_coverage = (send - sstart + 1) / slen
  ) %>%
  transmute(
    gene_name = qseqid,
    best_prior_hit = sseqid,
    pident, qlen, slen, aln_length = length,
    query_coverage, subject_coverage,
    reciprocal_best_hit = is_rbh,
    evalue, bitscore
  )

all_genes <- data.frame(gene_name = sub(" .*", "", names(readAAStringSet("venom_proteins.fasta"))))
rbh_table <- all_genes %>% left_join(rbh_table, by = "gene_name")

write.table(rbh_table, 'croVir4_croVir3_diamond_rbh_table.txt', col.names = T, row.names = F, quote = F, sep = '\t')

rbh_table <- rbh_table %>%
  mutate(novel_relative_to_prior_genome = case_when(
    is.na(reciprocal_best_hit) ~ "Yes",
    TRUE ~ "No"
  ))

evidence_matrix <- evidence_matrix %>%
  left_join(rbh_table, by = "gene_name")

colnames(evidence_matrix)
evidence_matrix <- evidence_matrix %>% 
  select(gene_name, seqnames, start, end, strand, n_exons, family, protein_length, 
         has_internal_stop, sp_prediction, sp_prob, has_signal_peptide, n_isoseq_transcripts, isoseq_supported, novel_relative_to_prior_genome.y)

# Intron support from RNA-seq splice junctions (SJ.tab.out)
introns <- exons %>%
  group_by(gene_name) %>%
  arrange(start) %>%
  mutate(
    intron_start = lag(end) + 1,
    intron_end = start - 1,
    seqnames = as.character(seqnames)
  ) %>%
  filter(!is.na(intron_start)) %>%
  select(gene_name, seqnames, intron_start, intron_end, strand) %>%
  ungroup()

introns_gr <- GRanges(introns$seqnames,
                      IRanges(introns$intron_start, introns$intron_end),
                      strand = introns$strand)

sj <- read.delim("CV1297_SJ.out.tab", header = FALSE,
                 col.names = c("seqnames", "intron_start", "intron_end", "strand_code",
                               "motif", "annotated", "n_unique_reads", "n_multi_reads", "max_overhang")) %>%
  mutate(strand = case_when(strand_code == 1 ~ "+", strand_code == 2 ~ "-", TRUE ~ "*")) %>%
  filter(n_unique_reads >= 1)

sj_gr <- GRanges(sj$seqnames, IRanges(sj$intron_start, sj$intron_end), strand = sj$strand)

hits <- findOverlaps(introns_gr, sj_gr, type = "equal", ignore.strand = TRUE)

introns$supported <- FALSE
introns$supported[queryHits(hits)] <- TRUE

intron_support_summary <- introns %>%
  group_by(gene_name) %>%
  summarise(pct_introns_supported = 100 * mean(supported), .groups = "drop")

evidence_matrix <- evidence_matrix %>%
  left_join(intron_support_summary, by = "gene_name")

evidence_matrix <- evidence_matrix %>% 
  select(-n_introns, -n_introns_supported)

table(evidence_matrix$pct_introns_supported)

write.table(evidence_matrix, 'croVir4_venoms_evidence_table.txt', col.names = T, row.names = F, quote = F, sep = '\t')
