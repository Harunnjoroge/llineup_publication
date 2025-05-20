# Adapted with simplifications from the script GAARD_work/randomisations/Fst/deprecated_scripts/fst_randomisations.r

# In this script, we load genotype data that is downloaded from Malariagen using the [./genomewide_association/glmm/extract whole genome genotype.ipynb ]
# We load a csv file containing permutated control_phase data the output of the [./genomewide_association/fst/llineup_random_permuation.r]
# Fst between the different cohorts are calculated and saved as rds. These are visualised using the [./genomewide_association/fst/calculate_Fst_window_p_value.r]

##############################################################################################################################################################################
library(data.table)
library(RColorBrewer)
library(magrittr)
library(stringr)
library(reticulate)
library(tidyverse)

chrom.order <- c('2L', '2R', '3L', '3R', 'X')

#Load snps data, split snp.id two chrom and pos and subset
setwd("/home/harunnn/lstm_scratch/network_scratch/llineup/llineup-genomics/llineup_publication/Data")
#load a subset of pops you want to analse differentiation 
gt_east_nonpbo <- './gt_east_nonpbo.csv'
gt_east_pbo <- './gt_east_pbo.csv'
gt_west_nonpbo <-'./gt_west_nonpbo.csv'
gt_west_pbo <-'./gt_west_pbo.csv'

#function to load data and split snps id 
load_data <- function(gt){
  df <- fread(gt,na.strings = '')
  df[,V1 :=NULL]
  colnames(df)[1] <-'snp_id'
  df[,c('Chrom', 'Pos'):=tstrsplit(snp_id, ':', fixed = TRUE, keep = 1:2)]
  df[,snp_id :=NULL]
  columns_to_move=c('Chrom', 'Pos')
  setcolorder(df, c(columns_to_move, setdiff(names(df), columns_to_move)))
  df[, Pos := as.integer(Pos)]
  setkey(df, Chrom, Pos)
  df[chrom.order]
}
gt <- list(east_nonpbo=gt_east_nonpbo, east_pbo=gt_east_pbo, west_nonpbo=gt_west_nonpbo, west_pbo = gt_west_pbo)
snp.tables <- lapply(gt,load_data)
populations <- names(snp.tables) 


#perform Fst analysis comparing both east and west uganda populations before net distribution and after

num.randomisations <-1000

# Load a conda env to run the reticulate code (which is inside the windowed.fst function)
use_condaenv('llineup')
allel <- import('allel')

# Set the window size for Fst analysis
window.size <- 1000


# Load the phenotype randomisations
setwd("./llineup_publication/Data")
phenotype.filename <- './net_round_randomisations_pop.csv'  # generate using the llineup_random_permutation.r script
phenotype.table <- fread(phenotype.filename, key = 'sample_id')
randomisations <- colnames(phenotype.table)[grepl('^r\\d+$', colnames(phenotype.table))]

samples.by.pop <- with(phenotype.table, split(sample_id, population))

#snp_table names
calculate.window.pos <- function(pop){
  # Get the number of windows that will be found on each chromosome
  chrom.snp.counts <- tapply(snp.tables[[pop]]$Chrom, snp.tables[[pop]]$Chrom, length)
  window.num <- ceiling(chrom.snp.counts / window.size)
  window.end <- cumsum(window.num)
  window.start <- window.end - window.num + 1
  snp.tables[[pop]][, window := rep(..window.start[.BY$Chrom]:..window.end[.BY$Chrom], each = ..window.size, length.out = ..chrom.snp.counts[.BY$Chrom]), by = 'Chrom']
  window.pos <- snp.tables[[pop]][, .(window.chrom = unique(Chrom), window.pos = floor(median(Pos))), by = window][, window := NULL]
  window.pos
}

# Write a function to pull out a phenotype vector 
get.phenotype.iteration <- function(pop, iteration){
  cat('\t', iteration, '\n', sep = '')
  with(phenotype.table[population == pop, c('sample_id', get('iteration'))], setNames(get(iteration), sample_id))
}

window.pos <- lapply(setNames(nm = populations), calculate.window.pos)
window.num.table <- sapply(window.pos, function(x) table(x$window.chrom))

# A function to turn a table of genotypes into a table of allele counts. 
get.allele.counts <- function(genotypes){
  mut.counts <- as.integer(rowSums(genotypes))
  wt.counts <- as.integer(ncol(genotypes) * 2 - mut.counts)
  data.table('wt'= wt.counts, 'mut'= mut.counts)
}

# Calculate the total allele count in each population, this will speed things up later

allele.counts.total <- lapply(setNames(nm = populations), function(pop) snp.tables[[pop]][, get.allele.counts(.SD), .SDcols = samples.by.pop[[pop]]])

# A function to calculate windowed Fst. Same as for the main data, except an extra element of parallelisation.
windowed.fst <- function(pop, phenotypes, window.size){
  
  samples.pre <- names(phenotypes)[phenotypes == 'pre']
  
  cat('\t\tCalculating allele counts.\n')
  allele.counts.pre <- snp.tables[[pop]][, .(Chrom, get.allele.counts(.SD)), .SDcols = samples.pre]
  # For the post, we can just subtract the pre from the total, saving us a bit of time. 
  allele.counts.post <- data.table(Chrom = snp.tables[[pop]]$Chrom, wt = allele.counts.total[[pop]]$wt - allele.counts.pre$wt, mut = allele.counts.total[[pop]]$mut - allele.counts.pre$mut)
  
  # Calculate moving Fst. 
  cat('\t\tCalculating Fst.\n')
  moving.Fst <- lapply(setNames(nm = chrom.order), function(chrom) matrix(allel$moving_patterson_fst(allele.counts.pre[Chrom == chrom, .(wt, mut)], allele.counts.post[Chrom == chrom, .(wt, mut)], as.integer(window.size)))) %>%
    # Pad with NAs where the final window of a chrom had < 1000 SNPs
    {lapply(names(.), 
            function(chrom) c(.[[chrom]], rep(NA, window.num.table[chrom, pop] - length(.[[chrom]]))))} 
  
  windowed.data <- data.table(moving.Fst = unlist(moving.Fst))
  windowed.data
}



windowed.fst.wrapper <- function(pop, randomisations, window.size){
  cat('\nRunning', length(randomisations), 'fst calculations for', pop, '\n')
  lapply(setNames(nm = c('control_phase', randomisations)), function(r) windowed.fst(pop, get.phenotype.iteration(pop, r), window.size))
}

randomised.windowed.data <- lapply(setNames(nm = populations), windowed.fst.wrapper, randomisations = randomisations[1:num.randomisations], window.size = window.size)

for (pop in names(randomised.windowed.data)){
  full.table <- do.call(cbind, randomised.windowed.data[[pop]])
  fst.table <- cbind(window.pos[[pop]], full.table[, grepl('Fst$', colnames(full.table)), with = F])
  saveRDS(fst.table, paste(pop, '_randomised_Fst.rds', sep = ''))
}

