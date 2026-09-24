paf_path    <- "old_to_new_full_genome.paf"
min_aln_len <- 2000
min_mapq    <- 20
gap_frac    <- 0.02  

chr_filter_regex <- "mi"   

paf_cols <- c("qname", "qlen", "qstart", "qend", "strand",
              "tname", "tlen", "tstart", "tend",
              "nmatch", "alnlen", "mapq")

paf_raw <- read_tsv(paf_path, col_names = FALSE, col_types = cols(.default = "c"),
                    progress = FALSE)
paf <- paf_raw[, 1:12]
names(paf) <- paf_cols
paf <- paf %>%
  mutate(across(c(qlen, qstart, qend, tlen, tstart, tend, nmatch, alnlen, mapq),
                as.numeric))

paf <- paf %>% filter(grepl(chr_filter_regex, qname), grepl(chr_filter_regex, tname))
if (nrow(paf) == 0) stop("No alignments have both qname and tname matching '", chr_filter_regex, "'.")

paf_f <- paf %>% filter(alnlen >= min_aln_len, mapq >= min_mapq)

## ---------------------------- Chromosome length tables -----------------------
old_lens <- paf %>% distinct(qname, qlen) %>% rename(chr = qname, len = qlen)
new_lens <- paf %>% distinct(tname, tlen) %>% rename(chr = tname, len = tlen)

## ---------------------------- Chromosome ordering -----------------------------
chrom_sort_key <- function(names) {
  base <- gsub("^scaffold-?", "", names, ignore.case = TRUE)
  is_mi <- grepl("^mi[0-9]+$", base, ignore.case = TRUE)
  num <- suppressWarnings(as.numeric(gsub("^mi", "", base, ignore.case = TRUE)))
  num[!is_mi] <- NA
  data.frame(type_rank = ifelse(is_mi, 1, 2), num = num)
}

order_chrom_table <- function(lens_df) {
  key <- chrom_sort_key(as.character(lens_df$chr))
  lens_df$type_rank <- key$type_rank
  lens_df$num <- key$num
  lens_df %>% arrange(type_rank, num, desc(len))
}

old_lens <- order_chrom_table(old_lens)
new_lens <- order_chrom_table(new_lens)

old_lens$chr <- factor(old_lens$chr, levels = old_lens$chr)
new_lens$chr <- factor(new_lens$chr, levels = new_lens$chr)

gap <- gap_frac * sum(old_lens$len)

make_offsets <- function(lens_df) {
  offsets <- numeric(nrow(lens_df))
  cum <- 0
  for (i in seq_len(nrow(lens_df))) {
    offsets[i] <- cum
    cum <- cum + lens_df$len[i] + gap
  }
  data.frame(chr = lens_df$chr, offset = offsets, len = lens_df$len)
}

old_off <- make_offsets(old_lens)
new_off <- make_offsets(new_lens)

paf_f <- paf_f %>%
  left_join(old_off %>% rename(qname = chr, q_offset = offset), by = "qname") %>%
  left_join(new_off %>% rename(tname = chr, t_offset = offset), by = "tname") %>%
  mutate(
    x1_top = q_offset + qstart,
    x2_top = q_offset + qend,
    x1_bot = ifelse(strand == "+", t_offset + tstart, t_offset + tend),
    x2_bot = ifelse(strand == "+", t_offset + tend,   t_offset + tstart)
  )

y_top <- 1
y_bot <- 0
ribbons <- paf_f %>%
  mutate(id = row_number()) %>%
  {
    df <- .
    data.frame(
      id    = rep(df$id, each = 4),
      x     = as.vector(rbind(df$x1_top, df$x2_top, df$x2_bot, df$x1_bot)),
      y     = rep(c(y_top, y_top, y_bot, y_bot), times = nrow(df)),
      qname = rep(df$qname, each = 4)
    )
  }

old_rects <- old_off %>% mutate(xmin = offset, xmax = offset + len)
new_rects <- new_off %>% mutate(xmin = offset, xmax = offset + len)

n_old <- nrow(old_lens)
old_colors <- setNames(hue_pal()(n_old), as.character(old_lens$chr))

ggplot() +
  geom_polygon(data = ribbons, aes(x = x, y = y, group = id, fill = qname),
               alpha = 0.55, color = NA) +
  geom_segment(data = old_rects, aes(x = xmin, xend = xmax, y = 1.02, yend = 1.02),
               linewidth = 4, color = "grey20") +
  geom_segment(data = new_rects, aes(x = xmin, xend = xmax, y = -0.02, yend = -0.02),
               linewidth = 4, color = "grey20") +
  geom_text(data = old_rects, aes(x = (xmin + xmax) / 2, y = 1.10, label = chr),
            angle = 30, hjust = 0, size = 3.6) +
  geom_text(data = new_rects, aes(x = (xmin + xmax) / 2, y = -0.10, label = chr),
            angle = -30, hjust = 1, size = 3.6) +
  scale_fill_manual(values = old_colors, guide = "none") +
  scale_y_continuous(limits = c(-0.45, 1.45)) +
  labs(title = "Synteny (microchromosomes only): croVir3 (top) vs croVir4 (bottom)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())
