library(data.table)
library(magrittr)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)

# Turn off scientific notation
options(scipen=999)

remove.inf <- function(x){
	unlist(x)[abs(unlist(x)) != 'Inf']
}

pal <- c(crash = 'goldenrod3', nocrash = 'royalblue1')


wf.filenames <- list.files('wf_simulations', 'W_values.*\\.csv')
wf <- data.table(W = numeric(), n = numeric(), N = numeric(), scenario = character())

for (fn in wf.filenames){
	N <- str_extract(fn, '(?<=_N)\\d+')
	this.wf <- fread(paste('wf_simulations', fn, sep = '/'))[, -c('V1')]
	for (cln in colnames(this.wf)){
		cln.split <- strsplit(cln, '_')[[1]]
		this.dt <- data.table(W = this.wf[[cln]], 
		                      n = as.numeric(sub('n', '', cln.split[3])),
		                      N = as.numeric(N),
		                      scenario = cln.split[2])
		wf <- rbind(wf,this.dt)
	}
}

coal.filenames <- list.files('coalescent_simulations', 'W_values.*\\.csv')
coal <- data.table(W = numeric(), n = numeric(), N = numeric(), scenario = character())

for (fn in coal.filenames){
	N <- str_extract(fn, '(?<=_N)\\d+')
	this.coal <- fread(paste('coalescent_simulations', fn, sep = '/'))[, -c('V1')]
	for (cln in colnames(this.coal)){
		cln.split <- strsplit(cln, '_')[[1]]
		this.dt <- data.table(W = this.coal[[cln]], 
		                      n = as.numeric(sub('n', '', cln.split[3])),
		                      N = as.numeric(N),
		                      scenario = cln.split[2])
		coal <- rbind(coal,this.dt)
	}
}

plot.median.and.quants <- function(W.table, title = '', return.plot = T) {
	if (nrow(W.table) == 0)
		return(ggplot(W.table, aes(x = n, y = W)) + geom_blank() + theme(panel.background = element_rect(fill = 'white')))
	medians.and.quants <- W.table[, as.list(setNames(quantile(W, probs = c(0.025, 0.5, 0.975), na.rm = T), c('Q.025', 'W', 'Q.975'))),
                           by = c('n','scenario')
                      ]
	y.lims <- c(0,
	            max(remove.inf(medians.and.quants[, 3:ncol(medians.and.quants)]), na.rm = T))
	x.lims <- c(0,500)
	# Now draw the actual plot
	pp <- ggplot(medians.and.quants, aes(x = n, y = W, color = scenario, fill = scenario)) + 
	             scale_fill_manual(values = pal) +
	             scale_colour_manual(values = pal) +
                 labs(y = 'Ne',title=title) +
                 geom_polygon(data = rbind(medians.and.quants[, .(n, scenario, conf = Q.025)], 
                                           medians.and.quants[nrow(medians.and.quants):1, .(n, scenario, conf = Q.975)]
                                     ), 
                              aes(x = n, y = conf),
                              alpha = 0.3,
                              size = 0) +
                 geom_line(size = 1.2) +
	             geom_point(pch = 21, size = 1.8, color = 'grey60') +
	             annotate('point', 
	                      x = medians.and.quants[Q.975 == Inf & scenario == 'crash', n],
	                      y = y.lims[2]*1.08, 
	                      pch = 8, size = 1.5, color = pal['crash']) +
	             annotate('point', 
	                      x = medians.and.quants[Q.975 == Inf & scenario == 'nocrash', n],
	                      y = y.lims[2]*1.13, 
	                      pch = 8, size = 1.5, color = pal['nocrash']) +
	             coord_cartesian(ylim = y.lims, xlim = x.lims, clip = 'off') +
                 theme_classic() +
	             theme(legend.title = element_blank())+
                 theme(plot.title = element_text(hjust = 0.5, vjust = 4)) 
	
	if (return.plot){
		return(pp)
	}
	else{
		print(pp)
		return(medians.and.quants)
	}
}

N.values <- sort(unique(coal$N))
wf.combined.plots <- lapply(N.values, function(pop) plot.median.and.quants(wf[N == pop], paste('Ne =', pop)))
coal.combined.plots <- lapply(N.values, function(pop) plot.median.and.quants(coal[N == pop], paste('Ne =', pop)))
combined.plots <- c(wf.combined.plots, coal.combined.plots)
combined.legend <- gtable_filter(ggplotGrob(combined.plots[[1]]), "guide-box")

svg('Ne_estimations.svg', width = 8, height = 12)
grid.arrange(
	do.call(arrangeGrob, c(lapply(combined.plots, function(p) p + theme(legend.position='none') + labs(x = '', y = '')), 
	                       nrow = 4,
	                       as.table = F,
	                       top = list(textGrob('', gp = gpar(fontface = "bold", cex = 1.3))),
	                       left = list(textGrob('Estimated Ne', rot = 90, vjust = 1)),
	                       bottom = list(textGrob('Sample size', vjust = -1)))
	),
	combined.legend,
	widths=unit.c(unit(1, "npc") - combined.legend$width, combined.legend$width), 
	nrow = 1, padding = unit(0, 'line')
)
dev.off()


r2 <- fread('power_calculations.csv') %>%
      .[, power := successes/trials]

pal2 <- setNames(rgb(seq(0.1, 0.4, 0.1), seq(0.4, 1, 0.2), seq(0.1, 0.4, 0.1), alpha = 0.7),
                 c('100', '1000', '10000', '100000')
)

plot.r2 <- function(r2.table, title = '') {
	x.lims <- c(0,500)
	y.lims <- c(0,1) 
	pp <- ggplot(r2.table, aes(x = n, y = power, color = as.character(N), fill = as.character(N))) + 
                 labs(y = 'Proportion',title=title) +
	             scale_colour_manual(values = pal2) +
	             scale_fill_manual(values = pal2) +
                 geom_line(size = 1.2) +
	             geom_point(pch = 21, size = 1.8, color = 'grey60') +
	             coord_cartesian(xlim = x.lims, ylim = y.lims, clip = 'off') +
                 theme_classic() +
	             labs(fill = 'Starting Ne', color = 'Starting Ne') +
                 theme(plot.title = element_text(hjust = 0.5, vjust = 4)) 
	
	return(pp)
}

combined.plots.r2 <- list(plot.r2(r2[simulation_type == 'wf'], title = 'Wright-Fisher'), 
                          plot.r2(r2[simulation_type == 'coalescent'], title = 'Coalescent'))
combined.legend.r2 <- gtable_filter(ggplotGrob(combined.plots.r2[[1]]), "guide-box")

svg('Ne_comparison_power.svg', width = 8, height = 4)
grid.arrange(
	do.call(arrangeGrob, c(lapply(combined.plots.r2, function(p) p + theme(legend.position='none') + labs(x = '', y = '')), 
	                       nrow = 1,
	                       top = list(textGrob('', gp = gpar(fontface = "bold", cex = 1.3))),
	                       left = list(textGrob('Proportion', rot = 90, vjust = 1)),
	                       bottom = list(textGrob('Sample size', vjust = -1)))
	),
	combined.legend.r2,
	widths=unit.c(unit(1, "npc") - combined.legend.r2$width, combined.legend.r2$width), 
	nrow = 1
)
dev.off()


