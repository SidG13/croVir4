# Load libraries ----
library(tidyverse)
library(Biostrings)

# Load data and format data ----
# Load in the landscape file
landscape_file <- readLines("data/croVir4_genome_2025-10-05_v1.Full_Mask.landscape")

# Find the header line for the Kimura portion of the file
kimura_idx <- which(grepl("^Class\\s+Repeat\\s+absLen", landscape_file))
coverage_idx <- which(grepl("^Coverage for each repeat class", landscape_file))

# Slice only between header+2 and the line before the coverage section
kimura_lines <- landscape_file[(kimura_idx + 2):(coverage_idx - 1)]

# Drop empty lines
kimura_lines <- kimura_lines[nzchar(trimws(kimura_lines))]

# Parse using fixed-width or whitespace splitting and make data frame
kimura_df <- bind_rows(
  lapply(
    kimura_lines,
    function(line) {
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      if (length(parts) >= 5) {
        data.frame(
          class = parts[1],
          repeat_name = parts[2],
          abs_len = as.numeric(parts[3]),
          well_char_len = as.numeric(parts[4]),
          kimura_pct = ifelse(parts[5] == "----", NA_real_, as.numeric(parts[5])),
          stringsAsFactors = FALSE
        )
      }
    }
  )
) |>
  # Set a floor for the columns that should be minimum zero
  mutate(
    abs_len = pmax(abs_len, 0),
    well_char_len = pmax(well_char_len, 0),
    kimura_pct = pmax(kimura_pct, 0)
  )

# Histogram of absolute length values
hist(
  x = kimura_df$abs_len,
  xlab = "Length of repeat",
  main = "Distribution of repeat lengths",
  breaks = 100
)

# Check for duplicates repeat_names across classes
dup_df <- kimura_df |>
  group_by(repeat_name) |>
  filter(n() > 1) |>
  arrange(repeat_name) |>
  select(class, repeat_name, abs_len)


# Genome length ----
# Get the genome
genome <- readDNAStringSet("/Users/ballardk/Library/CloudStorage/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/_Crotalus_XSpecies_Integration/croVir4_genomics/croVir4_genome_2025-10-05_v1.fasta")

# Get the genome size
genome_size <- sum(width(genome))
print(genome_size)


# Calculate fraction ----
# Calculate the fraction for each class
kimura_df2 <- kimura_df |>
  # For the duplicates, I will take the longest duplicated sequence only
  group_by(repeat_name) |>
  slice_max(order_by = abs_len, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(class, repeat_name, abs_len, well_char_len) |>
  mutate(
    genome_percentage = (abs_len / genome_size) * 100
  ) |>
  # Calculate by class
  group_by(class) |>
  summarise(
    genome_percentage = sum(genome_percentage),
    abs_len = sum(abs_len),
    well_char_len = sum(well_char_len)
  ) |>
  arrange(desc(genome_percentage))

# Calculate total
totals_df <- kimura_df2 |>
  summarise(
    genome_percentage = sum(genome_percentage),
    abs_len = sum(abs_len),
    well_char_len = sum(well_char_len)
  ) |>
  mutate(class = "Total")

# Calculate unmasked percentage
unmasked_df <- totals_df |>
  mutate(
    genome_percentage = 100 - genome_percentage,
    abs_len = genome_size - abs_len,
    well_char_len = genome_size - well_char_len,
    class = "Unmasked"
  )

# Add the totals
kimura_df2 <- kimura_df2 |>
  bind_rows(totals_df) |>
  bind_rows(unmasked_df)

# Create a list of colors
colors <- c(
  "Structural_RNA" = "#99993D",
  "Satellite" = "#CCCC52",
  "Simple_repeat" = "#FFFF66",
  "SINE/U" = "#DFCDF9",
  "SINE/MIR" = "#D7B4F8",
  "SINE/Deu" = "#CE9BF7",
  "SINE/tRNA-RTE" = "#C481F5",
  "SINE/tRNA-Alu" = "#BF74F4",
  "SINE/tRNA" = "#B966F4",
  "SINE/Alu" = "#B358F3",
  "SINE/7SL" = "#AD49F2",
  "SINE/5S" = "#A637F1",
  "SINE" = "#9F1FF0",
  "Retroposon/SVA" = "#FF4D4D",
  "PLE" = "#ACD8E5",
  "LINE/CRE" = "#C1D9FF",
  "LINE/R2" = "#A3C5DE",
  "LINE/Dong-R4" = "#99B3D7",
  "LINE/Jockey-I" = "#8FA1CF",
  "LINE/R1" = "#848FC8",
  "LINE/L2" = "#625CB1",
  "LINE/Rex-Babar" = "#554BAA",
  "LINE/CR1" = "#483AA2",
  "LINE/RTE" = "#38299A",
  "LINE" = "#251792",
  "LINE/L1" = "#00008B",
  "LTR/ERVK" = "#90ED90",
  "LTR/ERV" = "#81DD80",
  "LTR/ERV1" = "#73CD70",
  "LTR" = "#65BD61",
  "LTR/ERVL" = "#57AE51",
  "LTR/Gypsy" = "#489E42",
  "LTR/Copia" = "#3A8F33",
  "LTR/Pao" = "#2A8024",
  "LTR/DIRS" = "#006400",
  "RC/Helitron" = "#FF00FF",
  "DNA/Dada" = "#FFCFBC",
  "DNA/Transib" = "#FF9972",
  "DNA/TcMar" = "#FF936C",
  "DNA/Sola" = "#FF8D65",
  "DNA/PiggyBac" = "#FF865E",
  "DNA/P" = "#FF7F57",
  "DNA/MULE" = "#FF7850",
  "DNA/Merlin" = "#FF7149",
  "DNA" = "#FF6A42",
  "DNA/Maverick" = "#FF623B",
  "DNA/Kolobok" = "#FF5A34",
  "DNA/hAT" = "#FF512D",
  "DNA/Harbinger" = "#FF4825",
  "DNA/Ginger" = "#FF3D1E",
  "DNA/Crypton" = "#FF3115",
  "DNA/CMC" = "#FF200B",
  "Other" = "#4D4D4D",
  "Unknown" = "#999999",
  "Unmasked" = "black"
)

# Make data for plot
kimura_df3 <- kimura_df2 |>
  filter(class != "Total") |>
  # Recode the class to be less granular
  mutate(
    class_grouped = case_when(
      # DNA transposons
      str_detect(class, "^DNA/hAT") ~ "DNA/hAT",
      str_detect(class, "^DNA/TcMar") | class == "Mariner_Tc1" ~ "DNA/TcMar",
      str_detect(class, "^DNA/CMC") ~ "DNA/CMC",
      str_detect(class, "^DNA/MULE") ~ "DNA/MULE",
      str_detect(class, "^DNA/Kolobok") ~ "DNA/Kolobok",
      str_detect(class, "^DNA/PIF-Harbinger") ~ "DNA/Harbinger",
      str_detect(class, "^DNA/Crypton") ~ "DNA/Crypton",
      str_detect(class, "^DNA/EnSpm") | class == "DNA/En-Spm" | class == "DNA/CMC-EnSpm" ~ "DNA/CMC",
      class == "DNA/MuDR" ~ "DNA/MULE",
      class %in% c("DNA/Maverick", "DNA/Polinton") ~ "DNA/Maverick",
      class == "DNA/Transib" ~ "DNA/Transib",
      class == "DNA/Ginger" ~ "DNA/Ginger",
      class == "DNA/Merlin" ~ "DNA/Merlin",
      class == "DNA/Sola" ~ "DNA/Sola",
      class == "DNA/PiggyBac" ~ "DNA/PiggyBac",
      class == "DNA/P" | class == "DNA/P-Fungi" ~ "DNA/P",
      class == "DNA/Dada" ~ "DNA/Dada",
      class %in% c(
        "DNA/DNA", "DNA/DNA-7", "DNA/DNA-TA",
        "DNA/Chapaev", "DNA/T2", "DNA/Tc1",
        "DNA/SPIN", "DNA/SETARIA", "DNA/MERMITE",
        "DNA/Helitron", "DNA/CASINA", "DNA/Academ-1",
        "DNA-hAT", "DNA"
      ) ~ "DNA",
      # LTR
      str_detect(class, "^LTR/Gypsy") | class == "LTR/ty3/Gypsy" ~ "LTR/Gypsy",
      str_detect(class, "^LTR/ERV1") ~ "LTR/ERV1",
      str_detect(class, "^LTR/ERV4") ~ "LTR/ERV",
      str_detect(class, "^LTR/ERVK") ~ "LTR/ERVK",
      str_detect(class, "^LTR/ERVL") ~ "LTR/ERVL",
      str_detect(class, "^LTR/ERV-Foamy") | class == "LTR/ERV" ~ "LTR/ERV",
      str_detect(class, "^LTR/Copia") | class == "LTR_Copia" ~ "LTR/Copia",
      class == "LTR/Pao" ~ "LTR/Pao",
      class == "LTR/DIRS" ~ "LTR/DIRS",
      class == "LTR/BEL" ~ "LTR/Pao",
      class %in% c("LTR", "LTR/LTR") ~ "LTR",
      # LINE
      str_detect(class, "^LINE/CR1") ~ "LINE/CR1",
      str_detect(class, "^LINE/RTE") ~ "LINE/RTE",
      class == "LINE/L1" | class == "LINE/L1-Tx1" ~ "LINE/L1",
      class == "LINE/L2" | class == "LINE/L3" ~ "LINE/L2",
      class == "LINE/Rex-Babar" | class == "LINE/Rex1" ~ "LINE/Rex-Babar",
      class %in% c("LINE/Dong-R4", "LINE_R4/Dong", "LINE/R4") ~ "LINE/Dong-R4",
      class %in% c("LINE/Jockey", "LINE/I-Jockey") ~ "LINE/Jockey-I",
      class == "LINE/R1" ~ "LINE/R1",
      class == "LINE/R2" ~ "LINE/R2",
      class == "LINE/CRE" ~ "LINE/CRE",
      class %in% c(
        "LINE/Penelope", "LINEPenelope",
        "LINE/Tad1", "LINE/Daphne", "LINE/Kiri"
      ) ~ "LINE",
      class %in% c("LINE", "nonLTR", "nonLTR/brige") ~ "LINE",
      # SINE
      str_detect(class, "^SINE/tRNA-RTE") ~ "SINE/tRNA-RTE",
      str_detect(class, "^SINE/tRNA-Alu") ~ "SINE/tRNA-Alu",
      str_detect(class, "^SINE/tRNA") ~ "SINE/tRNA",
      class == "SINE/MIR" ~ "SINE/MIR",
      class == "SINE/Alu" ~ "SINE/Alu",
      class == "SINE/7SL" ~ "SINE/7SL",
      str_detect(class, "^SINE/5S") ~ "SINE/5S",
      class == "SINE/U" | class == "SINE/U-L1" ~ "SINE/U",
      str_detect(class, "^SINE/Deu") ~ "SINE/Deu",
      class %in% c(
        "SINE", "SINE/B2", "SINE/B4", "SINE/ID",
        "SINE/L2", "SINE/CORE", "SINE/SINE",
        "SINE/rRNA", "SINE/Sauria", "SINE-Sauria",
        "SINE-rRNA", "SINE?", "SINE?/L2"
      ) ~ "SINE",
      # PLE / Penelope
      str_detect(class, "^PENELOPE") | class == "PLE" ~ "PLE",
      # Retroposon
      class == "Retroposon/SVA" ~ "Retroposon/SVA",
      # RC
      class == "RC/Helitron" ~ "RC/Helitron",
      # Structural RNA
      class %in% c("tRNA", "rRNA", "snRNA", "scRNA", "RNA") ~ "Structural_RNA",
      # Satellite / Simple
      class == "Satellite" | str_detect(class, "^Satellite/") ~ "Satellite",
      class %in% c("Simple_repeat", "Simple_Repeat") ~ "Simple_repeat",
      # Catch-all
      class %in% c(
        "Unknown", "Unmasked", "Total",
        "Interspersed_Repeat", "Retroposon",
        "Retroposon/L1-dep", "Retroposon/L1-derived",
        "Retroposon/RTE-derived", "Retroposon/sno", "Retroposon/sno-L1"
      ) ~ class,
      TRUE ~ "Other"
    )
  ) |>
  group_by(class_grouped) |>
  summarise(
    genome_percentage = sum(genome_percentage),
    abs_len = sum(abs_len),
    .groups = "drop"
  )


# Genome fraction pie chart
(fraction_plt <- kimura_df2 |>
  filter(class != "Total") |>
  ggplot(aes(x = "", y = genome_percentage, fill = class)) +
  geom_col(color = "black", linewidth = 0.1) +
  coord_polar(theta = "y", start = 0) +
  labs(
    title = "Genome Fraction"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )
)
