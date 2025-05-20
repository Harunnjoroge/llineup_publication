library(data.table)
library(stringr)
library(RColorBrewer)
library(magrittr)

setwd("./llineup_publication/Data")

# complete llineup_seq data

meta <- fread('./merged_llineup_metadata.csv',na.strings = '')
meta <- unique(meta, by='sample_id')
#Remove first degree related samples fro KING relatedness analysis

samples.remove <- c('VBS50531-6645STDY11194268','VBS50533-6645STDY11194270','VBS50528-6645STDY11194265')
meta<- meta[!sample_id %in% samples.remove]
# change llin.actual entries
meta$LLIN.actual <- gsub("Non-PBO", "NonPBO", meta$LLIN.actual)

# Have a column to indicate population (location/insecticide)
meta$population <- with(meta, paste(LLIN.actual, sep = '_')) #randomise by round at each intervention seperately
#meta$population <- with(meta, paste('uganda')) #randomise by round the whole population

# Reorder and remove columns
column.order <- c('sample_id', 'Location', 'control_phase','population')
meta <- meta[, ..column.order]

#select round 

meta <- meta[control_phase %in% c('pre','post'), ] # ensure in your metadata you have a column where you have grouped your intervention to pre, intermediate and post so that you can select pre and post
# A function that will return a list of length k, where each k is a random shuffling of the input vector
shuffle <- function(x, k){
  data.table(replicate(k, sample(x, length(x), replace = F)))
}
# Create 1000 or more randomisations of the phenotype labels, stratified by population, and add them to the control
# table
set.seed(42)
num.randomisations <- 1000
replicate.names <- paste('r', str_pad(1:num.randomisations, nchar(as.character(num.randomisations)), pad = 0), sep = '')

meta[, (replicate.names) := shuffle(control_phase, num.randomisations), by = population]

# Write the table to file
fwrite(meta, './net_round_randomisations_pop.csv', sep = '\t', quote = F)

