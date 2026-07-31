library(tidyverse)
library(circlize)
library(ComplexHeatmap)

setwd('/home/administrator/ExtraDrive3/Yannick/CroVir4')

sex_chroms = c('scaffoldZ', 'scaffoldW')
scaffsize = read.delim('./chrom.list', header = FALSE, col.names = c('Chrom', 'size')) %>%
	arrange(Chrom %in% sex_chroms, match(Chrom, sex_chroms)) %>%
	mutate(
		is_micro = grepl('mi', Chrom),
		color = if_else(
			is_micro,
			if_else(cumsum(is_micro)  %% 2 == 1, "#f6dada", "#fcf0f0"),
			if_else(cumsum(!is_micro) %% 2 == 1, "#eaeaea", "#cbcbcb")
		)
	) %>%
	select(-is_micro) %>%
	mutate(Chrom = factor(Chrom, Chrom)) %>%
	as.data.frame()
scaffsize$color[scaffsize$Chrom == 'scaffoldZ'] = '#bffcf4'
scaffsize$color[scaffsize$Chrom == 'scaffoldW'] = '#f8fbd2'

gc = read.delim('./circos_tracks/gc_content.txt', header = FALSE,
	col.names = c('Chrom', 'start', 'end', 'gc')) %>%
	filter(Chrom %in% levels(scaffsize$Chrom))

repeats = read.delim('./circos_tracks/repeat_content.txt', header = FALSE,
	col.names = c('Chrom', 'start', 'end', 'repeat_frac')) %>%
	filter(Chrom %in% levels(scaffsize$Chrom))

genes = read.delim('./circos_tracks/gene_density.txt', header = FALSE,
	col.names = c('Chrom', 'start', 'end', 'nb_genes')) %>%
	filter(Chrom %in% levels(scaffsize$Chrom))

venom = read.delim('./venom_genes.txt', header = FALSE,
	col.names = c('Chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes')) %>%
	filter(feature == 'transcript') %>%
	mutate(
		gene   = gsub('.*gene_name "([^"]+)".*', '\\1', attributes),
		family = gsub('gene_name |"|Venom_|-.*', '', gene)
	) %>%
	select(Chrom, start, end, gene, family)

venom_family_pos = venom %>%
	arrange(Chrom, family, start) %>%
	group_by(Chrom, family) %>%
	mutate(group_id = cumsum(coalesce(start - lag(end), 0L) > 5e6)) %>%
	group_by(Chrom, family, group_id) %>%
	summarise(start = min(start), end = max(end), .groups = 'drop') %>%
	filter(!grepl('ADAM28', family)) %>%
	arrange(Chrom, start)

families = sort(unique(venom_family_pos$family))
n1 = ceiling(length(families) / 2)
n2 = length(families) - n1
set.seed(123)
family_colors = setNames(
	sample(c(colorspace::qualitative_hcl(n1, c = 85, l = 65),
	         colorspace::qualitative_hcl(n2, c = 70, l = 38))),
	families
)

venom_family_pos = venom_family_pos %>%
	mutate(col = family_colors[family])

gc      = gc      %>% mutate(col = ifelse(gc          > quantile(gc,          0.75), 'orange', 'black'))
repeats = repeats %>% mutate(col = ifelse(repeat_frac > quantile(repeat_frac, 0.75), 'orange', 'black'))

window = 500000
n_chrom = nrow(scaffsize)
col_gene = colorRamp2(c(0, 15, 30), c("white", "red", "#8b0000"))

pdf('circos.pdf', 12, 10)
circos.clear()
circos.par('track.height' = 0.2,
	cell.padding = c(0.01, 0, 0.01, 0),
	start.degree = 90,
	gap.degree = c(rep(1, n_chrom - 1), 20))
circos.initialize(scaffsize$Chrom, xlim = cbind(rep(1, n_chrom), scaffsize$size))

# Plot GC content
circos.track(ylim = range(gc$gc, na.rm = TRUE), track.height = 0.15, bg.border = NA, cell.padding = c(0,0,0,0),
	panel.fun = function(x, y) {

		c = scaffsize %>% filter(Chrom == CELL_META$sector.index)
		circos.rect(xleft = 0, xright = CELL_META$xlim[2],
			ybottom = CELL_META$ylim[1], ytop = CELL_META$ylim[2], col = c$color,
			border = NA)

		d = gc %>% filter(Chrom == CELL_META$sector.index) %>%
			mutate(run_id = cumsum(col != lag(col, default = first(col))))
		runs = split(seq_len(nrow(d)), d$run_id)
		for (col_target in c('black', 'orange')) {
			invisible(lapply(runs, function(idx) {
				if (d$col[idx[1]] != col_target) return(NULL)
				i = unique(c(max(1, min(idx) - 1), idx, min(nrow(d), max(idx) + 1)))
				circos.lines(x = d$start[i], y = d$gc[i], col = col_target)
			}))
		}

		circos.genomicAxis(h = 'top', labels.cex = par("cex"), tickLabelsStartFromZero = FALSE)

		circos.segments(
			x0 = CELL_META$xlim[2], y0 = CELL_META$ylim[2],
			x1 = CELL_META$xlim[2], y1 = CELL_META$ylim[2] + mm_y(5),
			lwd = 1)
		circos.text(CELL_META$xlim[2],
			CELL_META$cell.ylim[2] + mm_y(7),
			sprintf('%.0f Mb', CELL_META$xlim[2] / 1e6),
			facing = 'clockwise',
			niceFacing = TRUE,
			cex = 1,
			adj = c(0, 0.5))
})

# Plot repeat content
circos.track(ylim = range(repeats$repeat_frac, na.rm = TRUE), track.height = 0.15, bg.border = NA, cell.padding = c(0,0,0,0),
	panel.fun = function(x, y) {

		c = scaffsize %>% filter(Chrom == CELL_META$sector.index)
		circos.rect(xleft = 0, xright = CELL_META$xlim[2],
			ybottom = CELL_META$ylim[1], ytop = CELL_META$ylim[2], col = c$color,
			border = NA)

		d = repeats %>% filter(Chrom == CELL_META$sector.index) %>%
			mutate(run_id = cumsum(col != lag(col, default = first(col))))
		runs = split(seq_len(nrow(d)), d$run_id)
		for (col_target in c('black', 'orange')) {
			invisible(lapply(runs, function(idx) {
				if (d$col[idx[1]] != col_target) return(NULL)
				i = unique(c(max(1, min(idx) - 1), idx, min(nrow(d), max(idx) + 1)))
				circos.lines(x = d$start[i], y = d$repeat_frac[i], col = col_target)
			}))
		}
})

# Plot gene density
circos.track(ylim = c(0, 1), track.height = 0.05, cell.padding = c(0,0,0,0), bg.border = NA,
	panel.fun = function(x, y) {

		d = genes %>% filter(Chrom == CELL_META$sector.index) %>%
			mutate(col = col_gene(nb_genes))
		circos.rect(xleft = d$start, xright = d$start + window, ybottom = 0, ytop = 1, col = d$col, border = NA)
})

# Plot venom gene family positions
circos.track(ylim = c(0, 1), track.height = 0.06, cell.padding = c(0,0,0,0), bg.border = NA,
	panel.fun = function(x, y) {
		d = venom_family_pos %>% filter(Chrom == CELL_META$sector.index)
		if (nrow(d) == 0) return()
		circos.points(x = d$start, y = rep(CELL_META$ylim[2], nrow(d)), pch = 16, cex = 0.8, col = add_transparency(d$col, 0.35))
		circos.text(x = d$start, y = rep(CELL_META$ylim[2], nrow(d)),
			labels = d$family,
			facing = 'downward',
			cex = 0.7,
			adj = c(0, 1))
})

lgd_genes = Legend(col_fun = col_gene, title = "Gene density\n(genes / 500kb)",
	direction = "horizontal", title_position = "topcenter")
draw(lgd_genes, x = unit(0.82, "npc"), y = unit(0.08, "npc"))

dev.off()
