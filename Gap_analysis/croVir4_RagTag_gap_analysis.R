library(dplyr)
setwd('~/Documents/')

# Round 1 gaps (viridis pass)
r1_agp <- read.delim("ragtag.scaffold_r1.agp", header = FALSE, comment.char = "#",
                     col.names = c("object","object_start","object_end","part_num","component_type",
                                   "component_id","gap_length","gap_type","linkage","linkage_evidence"))
r1_gaps <- r1_agp %>%
  filter(component_type %in% c("N", "U")) %>%
  transmute(chrom = object, start = object_start, end = object_end,
            length = object_end - object_start + 1,
            round = "Round_1_viridis")

# Round 2 gaps (adamanteus pass)
r2_agp <- read.delim("ragtag.scaffold_r2.agp", header = FALSE, comment.char = "#",
                     col.names = c("object","object_start","object_end","part_num","component_type",
                                   "component_id","gap_length","gap_type","linkage","linkage_evidence"))
r2_gaps <- r2_agp %>%
  filter(component_type %in% c("N", "U")) %>%
  transmute(chrom = object, start = object_start, end = object_end,
            length = object_end - object_start + 1,
            round = "Round_2_adamanteus")

all_gaps <- bind_rows(r1_gaps, r2_gaps)

nrow(r1_gaps)
nrow(r2_gaps)  
nrow(all_gaps) 

all_gaps %>% count(chrom, round) %>% arrange(desc(n))

# --- NEW: map Round 1 objects to whichever Round 2 object absorbed them ---
r2_components <- r2_agp %>%
  filter(component_type == "W") %>%
  select(round2_object = object, component_id)

r1_gaps_mapped <- r1_gaps %>%
  left_join(r2_components, by = c("chrom" = "component_id")) %>%
  mutate(final_chrom = coalesce(round2_object, chrom))

r2_gaps_mapped <- r2_gaps %>%
  mutate(final_chrom = chrom)

combined_gaps <- bind_rows(
  r1_gaps_mapped %>% select(final_chrom, round, length),
  r2_gaps_mapped %>% select(final_chrom, round, length)
)

# --- NEW: rename Round 2 CM accessions to final published chromosome names ---
renaming_table <- c(
  "CM077918.1_RagTag" = "scaffoldma1", "CM077919.1_RagTag" = "scaffoldma2",
  "CM077920.1_RagTag" = "scaffoldma3", "CM077936.1_RagTag" = "scaffoldZ",
  "CM077921.1_RagTag" = "scaffoldma4", "CM077922.1_RagTag" = "scaffoldma5",
  "CM077923.1_RagTag" = "scaffoldma6", "CM077924.1_RagTag" = "scaffoldma7",
  "CM077930.1_RagTag" = "scaffoldmi1", "CM077926.1_RagTag" = "scaffoldmi2",
  "CM077925.1_RagTag" = "scaffoldmi3", "CM077928.1_RagTag" = "scaffoldmi4",
  "CM077935.1_RagTag" = "scaffoldW",   "CM077927.1_RagTag" = "scaffoldmi5",
  "CM077931.1_RagTag" = "scaffoldmi6", "CM077929.1_RagTag" = "scaffoldmi7",
  "CM077934.1_RagTag" = "scaffoldmi8", "CM077933.1_RagTag" = "scaffoldmi9",
  "CM077932.1_RagTag" = "scaffoldmi10"
)

combined_gaps <- combined_gaps %>%
  mutate(final_chrom_named = coalesce(renaming_table[final_chrom], final_chrom))

combined_gaps %>%
  count(final_chrom_named, round) %>%
  arrange(desc(n))

combined_gaps %>%
  group_by(final_chrom_named) %>%
  summarise(total_gaps = n(),
            round1_gaps = sum(round == "Round_1_viridis"),
            round2_gaps = sum(round == "Round_2_adamanteus"),
            .groups = "drop") %>%
  arrange(desc(total_gaps))
