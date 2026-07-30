library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

setwd('/Users/sidgopalan/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/_Crotalus_XSpecies_Integration/croVir4_genomics')

gtf <- rtracklayer::import('croVir4_annotation_2026-04-30_v2.gtf')
gtf <- as.data.frame(gtf)

# --- 1️⃣ Compute CDS lengths per transcript ---
cds_lengths <- gtf %>%
  filter(type == "CDS") %>%
  group_by(gene_name, transcript_id) %>%
  summarise(cds_length = sum(width), .groups = "drop")

# --- 2️⃣ Compute exon lengths per transcript ---
exon_lengths <- gtf %>%
  filter(type == "exon") %>%
  group_by(gene_name, transcript_id) %>%
  summarise(exon_length = sum(width), .groups = "drop")

# --- 3️⃣ Combine lengths ---
tx_lengths <- full_join(exon_lengths, cds_lengths,
                        by = c("gene_name", "transcript_id"))

# Replace NA CDS with 0 for easier logic
tx_lengths <- tx_lengths %>%
  mutate(cds_length = ifelse(is.na(cds_length), 0, cds_length))

# --- 4️⃣ Apply AGAT-like selection logic ---
longest_tx <- tx_lengths %>%
  group_by(gene_name) %>%
  group_modify(~{
    df <- .x
    
    if (any(df$cds_length > 0)) {
      # At least one isoform has CDS → choose longest CDS
      df %>% slice_max(order_by = cds_length,
                       n = 1,
                       with_ties = FALSE)
    } else {
      # No CDS in any isoform → choose longest exon length
      df %>% slice_max(order_by = exon_length,
                       n = 1,
                       with_ties = FALSE)
    }
  }) %>%
  ungroup()

# --- 5️⃣ Subset original GTF ---
gtf_longest <- gtf %>%
  filter(transcript_id %in% longest_tx$transcript_id)

gtf_longest %>%
  filter(type == "transcript") %>%
  count(gene_name) %>%
  summary()

gtf_longest <- gtf_longest %>% 
  mutate(gene_id = gene_name) %>% 
  mutate(transcript_id = gene_name)

gtf_longest <- gtf_longest %>% 
  select(seqnames, source, type, start, end, score, strand, phase, gene_name, gene_id, transcript_id) 

gr <- makeGRangesFromDataFrame(gtf_longest, keep.extra.columns = T)

rtracklayer::export(gr, "croVir4_annotation_2026-04-30_v2_KEEPLONGESTISOFORM-1.gtf") # Remove the ## lines written by rtracklayer after
