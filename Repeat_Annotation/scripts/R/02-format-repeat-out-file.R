# Load libraries ----
library(tidyverse)
library(data.table)
library(R.utils)
library(tidytable)
library(GenomicRanges)
library(Biostrings)
library(grid)
library(gtable)
library(ggrepel)

# Crotalus landscape data ----
## Load croVir3 data and format ----
out_file <- "data/croVir4_genome_2025-10-05_v1.Full_Mask.out.gz"
file.exists(out_file)

# Read in the .out file from the final step of RepeatMasker
out_df <- fread(
  out_file,
  skip = 3,
  fill = TRUE,
  header = FALSE
)

# Column names
colnames(out_df) <- c(
  "sw_score", "percent_divergence", "percent_deletion", "percent_insertions",
  "chrom", "chrom_start", "chrom_end", "chrom_left", 
  "strand", "repeat_name", "class", "repeat_start", "repeat_end", "repeat_left",
  "id", "empty"
)
glimpse(out_df)

# Remove the epmty column
out_df <- out_df |>
  select(-empty)

# Build the genomic ranges
gr <- GenomicRanges::GRanges(
  seqnames = out_df$chrom,
  ranges = IRanges(
    start = pmin(out_df$chrom_start, out_df$chrom_end),
    end = pmax(out_df$chrom_start, out_df$chrom_end)
  ),
  class = out_df$class
)

# Merge the overlaps in the ranges
total_gr <- reduce(gr, ignore.strand = TRUE)

# Get the genome
genome4 <- Biostrings::readDNAStringSet("data/croVir4_genome_2025-10-05_v1.fasta")

# Calculate the length of the genome
genome_size4 <- sum(width(genome4))
print(genome_size4)
# [1] 1585872272

# Get the masked base pairs
masked_bp <- sum(width(total_gr))
print(masked_bp)
# [1] 829030959

# Get the masked percents
masked_percent <- masked_bp / genome_size4 * 100
print(masked_percent)
# [1] 52.27602

# Print unique classes
print(distinct(out_df, class), n = 1000)

# Calculate classes
class_df <- out_df |>
  mutate(
    # Calculate sub-classes
    sub_class = case_when(
      # Class I: LTRs
      str_detect(class, regex("copia", ignore_case = TRUE)) ~ "Ty1/Copia",
      str_detect(class, regex("gypsy", ignore_case = TRUE)) ~ "Ty3/Gypsy",
      str_detect(class, regex("^LTR", ignore_case = TRUE)) ~ "LTR_unknown",

      # Non-LTR retrotransposons
      str_detect(class, regex("^LINE", ignore_case = TRUE)) ~ "LINE",
      str_detect(class, regex("nonLTR", ignore_case = TRUE)) ~ "LINE",
      str_detect(class, regex("Penelope", ignore_case = TRUE)) ~ "LINE",
      str_detect(class, regex("^SINE", ignore_case = TRUE)) ~ "SINE",

      # Class II: DNA transposons
      str_detect(class, regex("cacta", ignore_case = TRUE)) ~ "CACTA",
      str_detect(class, regex("MULE|MuDR", ignore_case = TRUE)) ~ "MULE",
      str_detect(class, regex("TcMar|Tc1|Mariner|Pogo", ignore_case = TRUE)) ~ "Tc1_Mariner",
      str_detect(class, regex("hat", ignore_case = TRUE)) ~ "hAT",
      str_detect(class, regex("Harbinger", ignore_case = TRUE)) ~ "PIF-Harbinger",
      str_detect(class, regex("polinton|maverick", ignore_case = TRUE)) ~ "Polinton",
      str_detect(class, regex("helitron", ignore_case = TRUE)) ~ "Helitron",

      # Other repeates
      str_detect(class, regex("Satellite", ignore_case = TRUE)) ~ "Satellite",
      str_detect(class, regex("Simple", ignore_case = TRUE)) ~ "Simple_repeat",
      str_detect(class, regex("Low_complexity", ignore_case = TRUE)) ~ "Low_complexity",
      str_detect(class, regex("rRNA|tRNA|snRNA|scRNA|srpRNA", ignore_case = TRUE)) ~ "RNA",
      str_detect(class, regex("Unknown", ignore_case = TRUE)) ~ "Unknown",
      TRUE ~ "Other"
    )
  )

# Get the ranges for each class
gr_list <- split(
  GRanges(
    seqnames = class_df$chrom,
    ranges = IRanges(
      start = pmin(class_df$chrom_start, class_df$chrom_end),
      end = pmax(class_df$chrom_start, class_df$chrom_end)
    )
  ),
  class_df$sub_class
)

# Summarise the classes by creating non-overlaping ranges
class_summary_df <- map_dfr(
  names(gr_list),
  function(cl) {
    reduced <- reduce(gr_list[[cl]])
    tibble(
      class = cl,
      bp = sum(width(reduced)),
      percent = sum(width(reduced)) / genome_size4 * 100
    )
  }
)

# Create a total percentage
totals_df <- tibble(
  class = "Total",
  bp = masked_bp,
  percent = masked_percent
)

# Create a unmasked percentage
unmasked_df <- totals_df |>
  mutate(
    class = "Unmasked",
    bp = genome_size4 - bp,
    percent = 100 - percent
  )

# Concatenate
class_summary_df <- class_summary_df |>
  bind_rows(totals_df) |>
  bind_rows(unmasked_df) |>
  mutate(
    # Calculate upper level classes
    top_class = case_when(
      # Class I: LTRs
      str_detect(class, regex("Ty1/Copia", ignore_case = TRUE)) ~ "Class I: LTR-RTs",
      str_detect(class, regex("Ty3/Gypsy", ignore_case = TRUE)) ~ "Class I: LTR-RTs",
      str_detect(class, regex("LTR_unknown", ignore_case = TRUE)) ~ "Class I: LTR-RTs",

      # Non-LTR retrotransposons
      str_detect(class, regex("LINE", ignore_case = TRUE)) ~ "Non-LTR-RTs",
      str_detect(class, regex("SINE", ignore_case = TRUE)) ~ "Non-LTR-RTs",

      # Class II: DNA transposons
      str_detect(class, regex("CACTA", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("MULE", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("Tc1_Mariner", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("hAT", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("Harbinger", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("polinton|maverick", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("helitron", ignore_case = TRUE)) ~ "Class II: DNA Transposon",

      # Other repeates
      str_detect(class, regex("Satellite", ignore_case = TRUE)) ~ "Other",
      str_detect(class, regex("Simple", ignore_case = TRUE)) ~ "Other",
      str_detect(class, regex("Low_complexity", ignore_case = TRUE)) ~ "Other",
      str_detect(class, regex("RNA", ignore_case = TRUE)) ~ "Other",
      str_detect(class, regex("Unknown", ignore_case = TRUE)) ~ "Other",
      str_detect(class, regex("Other", ignore_case = TRUE)) ~ "Other",
      TRUE ~ "Totals"
    )
  ) |>
  select(top_class, class, bp, percent) |>
  arrange(top_class)
print(class_summary_df, n = 100)
# > print(class_summary_df, n = 100)
# # A tidytable: 20 × 4
#    top_class                class                 bp percent
#    <chr>                    <chr>              <int>   <dbl>
#  1 Class I: LTR-RTs         LTR_unknown     55682523  3.51  
#  2 Class I: LTR-RTs         Ty1/Copia        8921790  0.563 
#  3 Class I: LTR-RTs         Ty3/Gypsy       81554787  5.14  
#  4 Class II: DNA Transposon CACTA             864599  0.0545
#  5 Class II: DNA Transposon Helitron         3828022  0.241 
#  6 Class II: DNA Transposon MULE             1345560  0.0848
#  7 Class II: DNA Transposon PIF-Harbinger    9819385  0.619 
#  8 Class II: DNA Transposon Polinton         1189722  0.0750
#  9 Class II: DNA Transposon Tc1_Mariner     56039822  3.53  
# 10 Class II: DNA Transposon hAT            119241441  7.52  
# 11 Non-LTR-RTs              LINE           399039746 25.2   
# 12 Non-LTR-RTs              SINE            52058413  3.28  
# 13 Other                    Low_complexity   8340281  0.526 
# 14 Other                    Other           31283140  1.97  
# 15 Other                    RNA              1190107  0.0750
# 16 Other                    Satellite        4647667  0.293 
# 17 Other                    Simple_repeat   59640654  3.76  
# 18 Other                    Unknown         19393619  1.22  
# 19 Totals                   Total          829030959 52.3   
# 20 Totals                   Unmasked       756841313 47.7   


## Load croVir3 data and format -----
out_file2 <- "data/Cvv_2017_genome_with_myos.Full_Mask.out.gz"
file.exists(out_file2)

# Read in the .out file from the final step of RepeatMasker
out_df2 <- fread(
  out_file2,
  skip = 3,
  fill = TRUE,
  header = FALSE
)

# Column names
colnames(out_df2) <- c(
  "sw_score", "percent_divergence", "percent_deletion", "percent_insertions",
  "chrom", "chrom_start", "chrom_end", "chrom_left", 
  "strand", "repeat_name", "class", "repeat_start", "repeat_end", "repeat_left",
  "id", "empty"
)
glimpse(out_df2)

# Remove the epmty column
out_df2 <- out_df2 |>
  select(-empty)

# Build the genomic ranges
gr2 <- GenomicRanges::GRanges(
  seqnames = out_df2$chrom,
  ranges = IRanges(
    start = pmin(out_df2$chrom_start, out_df2$chrom_end),
    end = pmax(out_df2$chrom_start, out_df2$chrom_end)
  ),
  class = out_df2$class
)

# Merge the overlaps in the ranges
total_gr2 <- reduce(gr2, ignore.strand = TRUE)

# Get the genome
genome3 <- Biostrings::readDNAStringSet("data/Cvv_2017_genome_with_myos.fasta")

# Calculate the length of the genome
genome_size3 <- sum(width(genome3))
print(genome_size3)
# [1] 1340207157

# Get the masked base pairs
masked_bp2 <- sum(width(total_gr2))
print(masked_bp2)
# [1] 551648355

# Get the masked percents
masked_percent2 <- masked_bp2 / genome_size3 * 100
print(masked_percent2)
# [1] 41.16142

# Print unique classes
print(distinct(out_df2, class), n = 1000)

# Calculate classes
class_df2 <- out_df2 |>
  mutate(
    # Calculate sub-classes
    sub_class = case_when(
      # Class I: LTRs
      str_detect(class, regex("copia", ignore_case = TRUE)) ~ "Ty1/Copia",
      str_detect(class, regex("gypsy", ignore_case = TRUE)) ~ "Ty3/Gypsy",
      str_detect(class, regex("^LTR", ignore_case = TRUE)) ~ "LTR_unknown",

      # Non-LTR retrotransposons
      str_detect(class, regex("^LINE", ignore_case = TRUE)) ~ "LINE",
      str_detect(class, regex("nonLTR", ignore_case = TRUE)) ~ "LINE",
      str_detect(class, regex("Penelope", ignore_case = TRUE)) ~ "LINE",
      str_detect(class, regex("^SINE", ignore_case = TRUE)) ~ "SINE",

      # Class II: DNA transposons
      str_detect(class, regex("cacta", ignore_case = TRUE)) ~ "CACTA",
      str_detect(class, regex("MULE|MuDR", ignore_case = TRUE)) ~ "MULE",
      str_detect(class, regex("TcMar|Tc1|Mariner|Pogo", ignore_case = TRUE)) ~ "Tc1_Mariner",
      str_detect(class, regex("hat", ignore_case = TRUE)) ~ "hAT",
      str_detect(class, regex("Harbinger", ignore_case = TRUE)) ~ "PIF-Harbinger",
      str_detect(class, regex("polinton|maverick", ignore_case = TRUE)) ~ "Polinton",
      str_detect(class, regex("helitron", ignore_case = TRUE)) ~ "Helitron",

      # Other repeates
      str_detect(class, regex("Satellite", ignore_case = TRUE)) ~ "Satellite",
      str_detect(class, regex("Simple", ignore_case = TRUE)) ~ "Simple_repeat",
      str_detect(class, regex("Low_complexity", ignore_case = TRUE)) ~ "Low_complexity",
      str_detect(class, regex("rRNA|tRNA|snRNA|scRNA|srpRNA", ignore_case = TRUE)) ~ "RNA",
      str_detect(class, regex("Unknown", ignore_case = TRUE)) ~ "Unknown",
      TRUE ~ "Other"
    )
  )

# Get the ranges for each class
gr_list2 <- split(
  GRanges(
    seqnames = class_df2$chrom,
    ranges = IRanges(
      start = pmin(class_df2$chrom_start, class_df2$chrom_end),
      end = pmax(class_df2$chrom_start, class_df2$chrom_end)
    )
  ),
  class_df2$sub_class
)

# Summarise the classes by creating non-overlaping ranges
class_summary_df2 <- map_dfr(
  names(gr_list2),
  function(cl) {
    reduced <- reduce(gr_list2[[cl]])
    tibble(
      class = cl,
      bp = sum(width(reduced)),
      percent = sum(width(reduced)) / genome_size3 * 100
    )
  }
)

# Create a total percentage
totals_df2 <- tibble(
  class = "Total",
  bp = masked_bp2,
  percent = masked_percent2
)

# Create a unmasked percentage
unmasked_df2 <- totals_df2 |>
  mutate(
    class = "Unmasked",
    bp = genome_size3 - bp,
    percent = 100 - percent
  )

# Concatenate
class_summary_df2 <- class_summary_df2 |>
  bind_rows(totals_df2) |>
  bind_rows(unmasked_df2) |>
  mutate(
    # Calculate upper level classes
    top_class = case_when(
      # Class I: LTRs
      str_detect(class, regex("Ty1/Copia", ignore_case = TRUE)) ~ "Class I: LTR-RTs",
      str_detect(class, regex("Ty3/Gypsy", ignore_case = TRUE)) ~ "Class I: LTR-RTs",
      str_detect(class, regex("LTR_unknown", ignore_case = TRUE)) ~ "Class I: LTR-RTs",

      # Non-LTR retrotransposons
      str_detect(class, regex("LINE", ignore_case = TRUE)) ~ "Non-LTR-RTs",
      str_detect(class, regex("SINE", ignore_case = TRUE)) ~ "Non-LTR-RTs",

      # Class II: DNA transposons
      str_detect(class, regex("CACTA", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("MULE", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("Tc1_Mariner", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("hAT", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("Harbinger", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("polinton|maverick", ignore_case = TRUE)) ~ "Class II: DNA Transposon",
      str_detect(class, regex("helitron", ignore_case = TRUE)) ~ "Class II: DNA Transposon",

      # Other repeates
      str_detect(class, regex("Satellite", ignore_case = TRUE)) ~ "Other",
      str_detect(class, regex("Simple", ignore_case = TRUE)) ~ "Other",
      str_detect(class, regex("Low_complexity", ignore_case = TRUE)) ~ "Other",
      str_detect(class, regex("RNA", ignore_case = TRUE)) ~ "Other",
      str_detect(class, regex("Unknown", ignore_case = TRUE)) ~ "Other",
      str_detect(class, regex("Other", ignore_case = TRUE)) ~ "Other",
      TRUE ~ "Totals"
    )
  ) |>
  select(top_class, class, bp, percent) |>
  arrange(top_class)
print(class_summary_df2, n = 100)
# Old:
# > print(class_summary_df2, n = 100)
# # A tidytable: 20 × 4
#    top_class                class                 bp percent
#    <chr>                    <chr>              <int>   <dbl>
#  1 Class I: LTR-RTs         LTR_unknown     31410004  2.34  
#  2 Class I: LTR-RTs         Ty1/Copia        7470028  0.557 
#  3 Class I: LTR-RTs         Ty3/Gypsy       43241026  3.23  
#  4 Class II: DNA Transposon CACTA             916993  0.0684
#  5 Class II: DNA Transposon Helitron         2197848  0.164 
#  6 Class II: DNA Transposon MULE             1119364  0.0835
#  7 Class II: DNA Transposon PIF-Harbinger    9205506  0.687 
#  8 Class II: DNA Transposon Polinton          522387  0.0390
#  9 Class II: DNA Transposon Tc1_Mariner     56381062  4.21  
# 10 Class II: DNA Transposon hAT             78096837  5.83  
# 11 Non-LTR-RTs              LINE           250626205 18.7   
# 12 Non-LTR-RTs              SINE            34302167  2.56  
# 13 Other                    Low_complexity   4622659  0.345 
# 14 Other                    Other           20088135  1.50  
# 15 Other                    RNA              1029819  0.0768
# 16 Other                    Satellite        1925986  0.144 
# 17 Other                    Simple_repeat   29254517  2.18  
# 18 Other                    Unknown         13658474  1.02  
# 19 Totals                   Total          551648355 41.2   
# 20 Totals                   Unmasked       788558802 58.8   
# New:
#   > print(class_summary_df2, n = 100)
# # A tidytable: 20 × 4
#    top_class                class                 bp percent
#    <chr>                    <chr>              <int>   <dbl>
#  1 Class I: LTR-RTs         LTR_unknown     36506063  2.72  
#  2 Class I: LTR-RTs         Ty1/Copia        7533252  0.562 
#  3 Class I: LTR-RTs         Ty3/Gypsy       47408676  3.54  
#  4 Class II: DNA Transposon CACTA             870361  0.0649
#  5 Class II: DNA Transposon Helitron         2348765  0.175 
#  6 Class II: DNA Transposon MULE             1144726  0.0854
#  7 Class II: DNA Transposon PIF-Harbinger    7577159  0.565 
#  8 Class II: DNA Transposon Polinton          316030  0.0236
#  9 Class II: DNA Transposon Tc1_Mariner     53562915  4.00  
# 10 Class II: DNA Transposon hAT             83819762  6.25  
# 11 Non-LTR-RTs              LINE           258052498 19.3   
# 12 Non-LTR-RTs              SINE            39055217  2.91  
# 13 Other                    Low_complexity   4740157  0.354 
# 14 Other                    Other           22861473  1.71  
# 15 Other                    RNA              1034592  0.0772
# 16 Other                    Satellite        2878142  0.215 
# 17 Other                    Simple_repeat   31146245  2.32  
# 18 Other                    Unknown         15809984  1.18  
# 19 Totals                   Total          561996649 41.9   
# 20 Totals                   Unmasked       778210508 58.1  

## Combine the croVir3 and croVir4 data ----
# Concatenate the dataframes
combined_class_summary_df <- bind_rows(list("croVir4" = class_summary_df, "croVir3" = class_summary_df2), .id = "genome")
glimpse(combined_class_summary_df)

# Format the data and write the data to a file
combined_class_summary_df |>
  pivot_wider(
    names_from = genome,
    values_from = c(bp, percent)
  ) |>
  # Add a percent change row
  mutate(
    percent_change = percent_croVir4 - percent_croVir3
  ) |>
  write_csv(file = "data/repeat-percentages.csv")


## Plot stuff ----
# Order the classes
top_class_order <- c(
  "Class I: LTR-RTs",
  "Non-LTR-RTs",
  "Class II: DNA Transposon",
  "Other",
  "Totals"
)

# Make data frame for plots
plot_df <- combined_class_summary_df |>
  mutate(
    top_class = factor(top_class, levels = top_class_order)
  ) |>
  arrange(top_class, percent) |>
  mutate(
    class = factor(class, levels = unique(class))
  )

# Plot the classes
(classes_plt <- ggplot(plot_df, aes(x = class, y = percent, fill = top_class)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.75) +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1),
    limits = c(0, 65),
    expand = expansion()
  ) +
  scale_fill_brewer(
    palette = "Dark2"
  ) +
  geom_text(
    aes(label = paste0(round(percent, digits = 1), "%")),
    vjust = -0.75, size = 3
  ) +
  labs(
    title = "Repetitive sequence annotations",
    fill = "Type",
    x = NULL,
    y = "Genome %"
  ) +
  theme_minimal() +
  theme(
    aspect.ratio = 0.4,
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85),
    legend.background = element_rect(linewidth = 0.25)
  ) +
  facet_wrap(~genome, ncol = 1)
)
ggsave(
  plot = classes_plt,
  filename = "figures/repeat-percentages.pdf",
  width = 10, height = 8, create.dir = TRUE
)

# Plot the classes as dots with a line between them for each genome
(genome_comp_plt <- ggplot(plot_df,aes(x = genome, y = percent, color = genome, group = class)) +
  geom_path(color = "black") +
  geom_point() +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1),
    limits = c(0, 65)
  ) +
  ggrepel::geom_text_repel(
    aes(
      label = paste0(round(percent, digits = 1), "%"),
      vjust = if_else(genome == "croVir4", -0.75, 1.75)
    ),
    size = 3, show.legend = FALSE
  ) +
  coord_cartesian(clip = "off") +
  scale_color_brewer(
    palette = "Dark2",
    labels = c("croVir4" = "This study", "croVir3" = "Schield et al. 2019")
  ) +
  labs(
    title = "Comparison between genome versions",
    color = "Genome version",
    x = NULL,
    y = "Genome %"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(face = "bold"),
    # Remove spacing between the panels in the x-axis
    panel.spacing.x = unit(0, "pt"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(), angle = 80),
    # Remove x-axis ticks
    axis.ticks.x = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85),
    legend.background = element_rect(linewidth = 0.25)
  ) +
  facet_wrap(~class, nrow = 1, strip.position = "bottom")
)

# Put a border around all the panels together
# Convert to gtable and add a single border around all panels
genome_comp_plt <- ggplotGrob(genome_comp_plt)

# Find the extent of all panel cells
panels <- genome_comp_plt$layout[grepl("^panel", genome_comp_plt$layout$name), ]
top <- min(panels$t)
bottom <- max(panels$b)
left <- min(panels$l)
right <- max(panels$r)

# Add a rectangle grob spanning all panels
border_grob <- rectGrob(gp = gpar(col = "black", fill = NA, lwd = 0.25 * .pt))
genome_comp_plt <- gtable_add_grob(genome_comp_plt, border_grob, t = top, b = bottom, l = left, r = right, z = Inf, name = "panel-border")

# Re-draw the plot
grid.newpage()
grid.draw(genome_comp_plt)
ggsave(
  plot = genome_comp_plt,
  filename = "figures/repeat-percentages-comparison-between-genome-versions.svg",
  width = 10, height = 5, create.dir = TRUE
)
