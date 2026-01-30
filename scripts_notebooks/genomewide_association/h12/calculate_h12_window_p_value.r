#The script calculate delta H12 and visualise it.

library(stringi)
library(data.table)
library(magrittr)
library(future.apply)
plan(tweak(multisession, workers = 30))

# Load all of the H12 tables
setwd("~/lstm_scratch/network_scratch/llineup/llineup-genomics")
h12.filenames <- list.files('./h12/data', pattern = '\\.tsv$', full.names = TRUE)


study.pops <- setNames(nm = c('East_NonPBO', 'East_PBO', 'West_NonPBO',
                             'West_PBO'))

# Get the randomisation ids. 
randomisation.ids <- readLines(h12.filenames[1], n = 1) %>%
                     {strsplit(., '\t')[[1]]} %>%
                     grep('r\\d{5}', ., value = T)


# A function that looks for positive peaks by identifying windows more extreme than thrice the 
# 95% centile:
find.peaks <- function(values, centile = 0.95, multiplier = 3){
    center <- median(values)
	thresh1 <- center + (quantile(values, centile) - center) * multiplier
	thresh2 <- center + (quantile(values, 1-centile) - center) * multiplier
	(values > thresh1) | (values < thresh2)
}

# A function to load and combine all data for a given randomisation id
load.and.combine.data <- function(pop, find.peaks =T, calculate.P.values = T){
	
	pre.filenames <- grep(paste(pop, '\\.pre', sep = ''), h12.filenames, value = T) %>%
	                   setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.tsv)'))
	post.filenames <- grep(paste(pop, '\\.post', '\\.', sep = ''), h12.filenames, value = T) %>%
	                  setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.tsv)'))
	
	if (any(names(post.filenames) != names(pre.filenames)))
		stop(paste('Different chromosomes have been found for post and pre samplesets for randomisation', rand.id))
	
	pre.data <- names(pre.filenames) %>%
	              lapply(function(chrom){ 
	                  fread(pre.filenames[chrom], sep = '\t') %>%
	                  .[, chromosome := ..chrom] %>%
					  setnames('h12', 'H12') %>%
					  setcolorder(c('chromosome', 'startpoint', 'endpoint', 'midpoint'))
				  }) %>%
				  rbindlist()
	post.data <- names(post.filenames) %>%
	             lapply(function(chrom){ 
	                 fread(post.filenames[chrom], sep = '\t') %>%
	                 .[, chromosome := ..chrom] %>%
					 setnames('h12', 'H12') %>%
					 setcolorder(c('chromosome', 'startpoint', 'endpoint', 'midpoint'))
				 }) %>%
				 rbindlist()
	# Check that the midpoints are identical
	if (!identical(pre.data$midpoint, post.data$midpoint))
		stop('Window positions do no match between pre and post samples')
	diff.data <- cbind(post.data[, .(chromosome, startpoint, endpoint, midpoint)],
	                   post.data[, c('H12', ..randomisation.ids)] - pre.data[, c('H12', ..randomisation.ids)]
	)
	
	# Look for peaks in the diff data
	if (find.peaks){
		diff.data[, is.peak := find.peaks(H12, centile = 0.98)]
		setcolorder(diff.data, c('chromosome', 'startpoint', 'endpoint', 'midpoint', 'H12', 'is.peak'))
	}
	
	if (calculate.P.values){
		if(find.peaks){
			diff.data[is.peak == T, pval := apply(.SD - H12, 
			                                      1, 
			                                      # get TWO-TAILED p-value
			                                      function(x) 1-(abs(sum(x >= 0)/length(x) -0.5)*2)
			                                ), 
			                                .SDcols = randomisation.ids
			] 
		}
		else {
			diff.data[, pval := apply(.SD - H12, 
			                          1, 
			                          function(x) 1-(abs(sum(x >= 0)/length(x) -0.5)*2)
			                    ), 
			                    .SDcols = randomisation.ids
			] 
		}
	}
	
	output <- list(pre = pre.data, post = post.data, diff = diff.data)
	output
}

cat('Loading data\n')
h12.tables <- lapply(study.pops, function(pop) {cat(pop, '\n'); load.and.combine.data(pop)})

source("scripts/R_plotting.r")

chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

# Function to plot the H12 data
plot.h12.diff <- function(h12.table, 
                          filename = NULL, 
                          num.randomisations = NULL, 
                          p.thresh = 0.01, 
                          plot.title = '', 
                          gaps = 5000000, 
                          filter.name = 'is.peak'){
	# Create the plot
	if (!missing(filename)){
		file.width = 6.5
		file.height = 3.5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	# Get vectors of start and end points for each chromosome (ie: cumulative sizes + gaps)
	ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
	cs <- ce - chrom.sizes
	layout(matrix(c(rep(1,4),rep(2,1)), nrow = 5, ncol = 1))
	colours <- c(h12 = lighten.col('red', alpha = 0.8),
                 randomisations = lighten.col('grey80', alpha = 0.15))
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	num.randomisations <- ifelse(is.null(num.randomisations), 
	                             length(randomisation.ids),
	                             num.randomisations)
	# This odd way of getting a sequence is to make sure we get the right outcome if num.randomisations == 0
	r.columns <- randomisation.ids[seq(1, num.randomisations, length.out = num.randomisations)]
	h12.columns <- c('H12', r.columns)
	max.y <- max(c(max(h12.table[, ..h12.columns]), 0.05))
	min.y <- min(h12.table[, ..h12.columns])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	h12.table$genome.pos <- h12.table$midpoint + cs[as.character(h12.table$chromosome)]
	# Add randomised data
	sapply(r.columns, function(x) h12.table[, lines(genome.pos, get(x), col = colours['randomisations'], lwd = 0.4), by = chromosome])
	# Add true data
	h12.table[, lines(genome.pos, H12, col = colours['h12'], lwd = 1.2), by = chromosome]
	# Add y axis
	h12.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, h12.step.size))
	mtext('H12', 2, 2, cex = 0.8)#'H12 difference
	# Add peaks that pass filtering. Colour according to whether they have significant p-value
	if (length(p.thresh) == 1)
		p.thresh <- c(p.thresh, -1)
	p.colours <- c('orchid3', 'green', 'blue')
	filter.pass <- h12.table[[filter.name]]
	h12.table[filter.pass, points(genome.pos, H12, 
	                              pch = 21,
	                              bg = p.colours[(pval < p.thresh[1]) + (pval < p.thresh[2]) + 1], 
	                              col = colours['randomisations'], 
	                              cex = 1.1, lwd = .5)
	]
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(gaps = gaps, chrom.sizes, gene.cex = 0.7, point.cex = 1, chrom.offset = -1.2, chrom.cex = 1.2)
	
	if (!missing(filename))
		dev.off()
}
  
p.threshold = 0.01

for (pop in names(h12.tables))
	plot.h12.diff(h12.tables[[pop]]$diff, #pre, post or diff
	              filename = paste('./h12/data',pop, '_delta_h12.png'),
	              p.thresh = p.threshold,
	              plot.title = pop)




