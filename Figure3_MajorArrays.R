library(tidyverse)
library(rtracklayer)
library(RColorBrewer)
library(scales)
library(cowplot)
library(gggenes)
library(ggrepel)

setwd('/Volumes/SeagatePortableDrive/croVir4_Genome_Figs')


#### START Gene tracks for croVir3 and croVir4 ####
gtf_croVir3 <- rtracklayer::import('/Volumes/SeagatePortableDrive/snake_genomes/Crotalus_viridis/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod.gtf') %>% 
  as.data.frame()
gtf_croVir4 <- rtracklayer::import('/Volumes/SeagatePortableDrive/snake_genomes/Crotalus_viridis/croVir4_Annotation_2025/croVir4_annotation_2026-07-28_v1.gtf') %>% 
  as.data.frame()

### START croVir3 SVSP ###

SVSP.info_3 <- gtf_croVir3 %>% 
  filter(seqnames == 'scaffold-mi2' & start > 8503562 & end < 9023053) %>% 
  select(seqnames, start, end, strand, type, gene_id) %>% 
  filter(type == 'gene') %>% 
  filter(str_detect(gene_id, 'SVSP|810|808|833|823')) %>% # adding non venom bookends
  filter(str_detect(gene_id, 'fgenesh|maker')) %>% 
  mutate(strand = if_else(strand == '+', 1, -1)) %>% 
  mutate(seqnames = gsub('scaffold-','', seqnames)) %>% 
  mutate(gene_type = if_else(grepl('SVSP', gene_id),  'Ven', 'NonVen'))

SVSP.reg.start_3 <- 8503562 # SVSP
SVSP.reg.end_3 <- 9023053 # SVSP

SVSP.reg.length_3 <- paste(c(round((SVSP.reg.end_3-SVSP.reg.start_3)/1000,digits = 2),'kb'),collapse = ' ')

p_genes <- ggplot(SVSP.info_3, aes(xmin = start, xmax = end, y = 'gene_track', forward = strand, fill = gene_type)) +
  geom_segment(aes(x=SVSP.reg.start_3, xend=SVSP.reg.end_3, y='gene_track', yend='gene_track'), lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
  geom_text_repel(aes(x = (start + end)/2, label = sub(".*-", "", gene_id)), nudge_y = 0.2, size = 3, direction = "x") +
  ylab('') +
  xlab('') +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("Ven" = "darkgreen", "NonVen" = "grey60")) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank()) +
  scale_x_continuous(
    expand = c(0,0),
    labels = scales::label_number(scale = 1e-6),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(x = "Position on microchromosome 2 (Mb)") +
  theme_minimal()

### START croVir4 SVSP ###

SVSP.info_4 <- gtf_croVir4 %>% 
  filter(seqnames == 'scaffoldmi2' & start > 10986736 & end < 11645030) %>% 
  select(seqnames, start, end, strand, type, gene_name) %>% 
  filter(type == 'transcript') %>% 
  mutate(strand = if_else(strand == '+', 1, -1)) %>% 
  mutate(seqnames = gsub('scaffold-','', seqnames)) %>% 
  mutate(gene_type = if_else(grepl('SVSP', gene_name),  'Ven', 'NonVen'))

SVSP.reg.start_4 <- 10986736 # SVSP
SVSP.reg.end_4 <- 11645030 # SVSP

SVSP.reg.length_4 <- paste(c(round((SVSP.reg.end_4-SVSP.reg.start_4)/1000,digits = 2),'kb'),collapse = ' ')

p_genes2 <- ggplot(SVSP.info_4, aes(xmin = start, xmax = end, y = 'gene_track', forward = strand, fill = gene_type)) +
  geom_segment(aes(x=SVSP.reg.start_4, xend=SVSP.reg.end_4, y='gene_track', yend='gene_track'), lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
  geom_text_repel(aes(x = (start + end)/2, label = gsub('Venom_', '', gene_name)), nudge_y = 0.2, size = 3, direction = "x") +
  ylab('') +
  xlab('') +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("Ven" = "darkgreen", "NonVen" = "grey60")) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank()) +
  scale_x_continuous(
    expand = c(0,0),
    labels = scales::label_number(scale = 1e-6),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(x = "Position on microchromosome 2 (Mb)") +
  theme_minimal()

p1 <- plot_grid(p_genes, p_genes2, nrow = 2, axis = 'lr', align = 'hv') # (10.86 x 6.82 landscape)

### END SVSP ###

### START croVir3 SVMP ###

SVMP.info_3 <- gtf_croVir3 %>% 
  filter(seqnames == 'scaffold-mi1' & start > 13778830 & end < 14621367) %>% 
  select(seqnames, start, end, strand, type, gene_id) %>% 
  filter(type == 'gene') %>% 
  filter(str_detect(gene_id, 'SVMP|12843|12844|12864|12854|12856')) %>% # adding non venom bookends
  filter(str_detect(gene_id, 'fgenesh|maker')) %>% 
  mutate(strand = if_else(strand == '+', 1, -1)) %>% 
  mutate(seqnames = gsub('scaffold-','', seqnames)) %>% 
  mutate(gene_type = if_else(grepl('SVMP', gene_id),  'Ven', 'NonVen'))

SVMP.reg.start_3 <- 13778830 # SVSP
SVMP.reg.end_3 <- 14621367 # SVSP

SVMP.reg.length_3 <- paste(c(round((SVMP.reg.end_3-SVMP.reg.start_3)/1000,digits = 2),'kb'),collapse = ' ')

p_genes_3 <- ggplot(SVMP.info_3, aes(xmin = start, xmax = end, y = 'gene_track', forward = strand, fill = gene_type)) +
  geom_segment(aes(x=SVMP.reg.start_3, xend=SVMP.reg.end_3, y='gene_track', yend='gene_track'), lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
  geom_text_repel(aes(x = (start + end)/2, label = sub(".*-", "", gene_id)), nudge_y = 0.2, size = 3, direction = "x") +
  ylab('') +
  xlab('') +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("Ven" = "darkgreen", "NonVen" = "grey60")) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank()) +
  scale_x_continuous(
    expand = c(0,0),
    labels = scales::label_number(scale = 1e-6),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(x = "Position on microchromosome 1 (Mb)") +
  theme_minimal()

### START croVir4 SVMP ###

SVMP.info_4 <- gtf_croVir4 %>% 
  filter(seqnames == 'scaffoldmi3' & start > 14623826 & end < 15565670) %>% 
  select(seqnames, start, end, strand, type, gene_name) %>% 
  filter(type == 'transcript') %>% 
  mutate(strand = if_else(strand == '+', 1, -1)) %>% 
  mutate(seqnames = gsub('scaffold-','', seqnames)) %>% 
  mutate(gene_type = if_else(grepl('SVMP', gene_name),  'Ven', 'NonVen'))

SVMP.reg.start_4 <- 14623826 # SVMP
SVMP.reg.end_4 <- 15565670 # SVMP

SVMP.reg.length_4 <- paste(c(round((SVMP.reg.end_4-SVMP.reg.start_4)/1000,digits = 2),'kb'),collapse = ' ')

p_genes_4 <- ggplot(SVMP.info_4, aes(xmin = start, xmax = end, y = 'gene_track', forward = strand, fill = gene_type)) +
  geom_segment(aes(x=SVMP.reg.start_4, xend=SVMP.reg.end_4, y='gene_track', yend='gene_track'), lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
  geom_text_repel(aes(x = (start + end)/2, label = gsub('Venom_', '', gene_name)), nudge_y = 0.2, size = 3, direction = "x") +
  ylab('') +
  xlab('') +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("Ven" = "darkgreen", "NonVen" = "grey60")) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank()) +
  scale_x_continuous(
    expand = c(0,0),
    labels = scales::label_number(scale = 1e-6),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(x = "Position on microchromosome 3 (Mb)") +
  theme_minimal()

p2 <- plot_grid(p_genes_3, p_genes_4, nrow = 2, axis = 'lr', align = 'hv') # (10.86 x 6.82 landscape)

### START croVir3 PLA2 ###

PLA2.info_3 <- gtf_croVir3 %>% 
  filter(seqnames == 'scaffold-mi7' & start > 3011365 & end < 3056355) %>% 
  select(seqnames, start, end, strand, type, gene_id) %>% 
  filter(type == 'gene') %>% 
  filter(str_detect(gene_id, 'fgenesh|maker')) %>% 
  mutate(strand = if_else(strand == '+', 1, -1)) %>% 
  mutate(seqnames = gsub('scaffold-','', seqnames)) %>% 
  mutate(gene_type = if_else(grepl('PLA2', gene_id),  'Ven', 'NonVen'))

PLA2.reg.start_3 <- 3011365 # PLA2
PLA2.reg.end_3 <- 3056355 # PLA2

PLA2.reg.length_3 <- paste(c(round((PLA2.reg.end_3-PLA2.reg.start_3)/1000,digits = 2),'kb'),collapse = ' ')

p_genes_5 <- ggplot(PLA2.info_3, aes(xmin = start, xmax = end, y = 'gene_track', forward = strand, fill = gene_type)) +
  geom_segment(aes(x=PLA2.reg.start_3, xend=PLA2.reg.end_3, y='gene_track', yend='gene_track'), lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
  geom_text_repel(aes(x = (start + end)/2, label = sub(".*-", "", gene_id)), nudge_y = 0.2, size = 3, direction = "x") +
  ylab('') +
  xlab('') +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("Ven" = "darkgreen", "NonVen" = "grey60")) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank()) +
  scale_x_continuous(
    expand = c(0,0),
    labels = scales::label_number(scale = 1e-6),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(x = "Position on microchromosome 7 (Mb)") +
  theme_minimal()

### START croVir4 PLA2 ###

PLA2.info_4 <- gtf_croVir4 %>% 
  filter(seqnames == 'scaffoldmi6' & start > 9902610 & end < 9950209) %>% 
  select(seqnames, start, end, strand, type, gene_name) %>% 
  filter(type == 'transcript') %>% 
  mutate(strand = if_else(strand == '+', 1, -1)) %>% 
  mutate(seqnames = gsub('scaffold-','', seqnames)) %>% 
  mutate(gene_type = if_else(grepl('PLA2', gene_name),  'Ven', 'NonVen')) %>% 
  distinct()

PLA2.reg.start_4 <- 9902610
PLA2.reg.end_4 <- 9950209 

PLA2.reg.length_4 <- paste(c(round((PLA2.reg.end_4-PLA2.reg.start_4)/1000,digits = 2),'kb'),collapse = ' ')

p_genes_6 <- ggplot(PLA2.info_4, aes(xmin = start, xmax = end, y = 'gene_track', forward = strand, fill = gene_type)) +
  geom_segment(aes(x=PLA2.reg.start_4, xend=PLA2.reg.end_4, y='gene_track', yend='gene_track'), lwd=1,color='grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
  geom_text_repel(aes(x = (start + end)/2, label = gsub('Venom_', '', gene_name)), nudge_y = 0.2, size = 3, direction = "x") +
  ylab('') +
  xlab('') +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("Ven" = "darkgreen", "NonVen" = "grey60")) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank()) +
  scale_x_continuous(
    expand = c(0,0),
    labels = scales::label_number(scale = 1e-6),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(x = "Position on microchromosome 6 (Mb)") +
  theme_minimal()

p3 <- plot_grid(p_genes_5, p_genes_6, nrow = 2, axis = 'lr', align = 'hv') # (10.86 x 6.82 landscape)

plot_grid(p1, p2, p3, nrow = 3, axis = 'lr', align = 'hv')
