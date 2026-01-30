library(data.table)
library(stringr)
library(RColorBrewer)
library(magrittr)

setwd('~/lstm_scratch/network_scratch/llineup/llineup-genomics')

# complete llineup_seq data

meta <- fread('./data/merged_llineup_metadata.csv',na.strings = '')
meta <- unique(meta, by='sample_id')
#Remove first degree related samples fro KING relatedness analysis

samples.remove <- c('VBS50531-6645STDY11194268','VBS50533-6645STDY11194270','VBS50528-6645STDY11194265')
meta<- meta[!sample_id %in% samples.remove]
# # change llin.actual entries
# meta$LLIN.actual <- gsub("Non-PBO", "NonPBO", meta$LLIN.actual)
#select round 

meta <- meta[control_phase %in% c('pre','post'), ] # ensure in your metadata you have a column where you have grouped your intervention to pre, intermediate and post so that you can select pre and post


# Have a column to indicate population (here we use control_phase, but could be location or other grouping variable)
meta$population <- with(meta, paste(control_phase)) #randomise by intervention phase for fst between east and west(pre and post)
#for h12, you may want to randomise by location only so that you keep the pre and post structure within each location(able to detect changes in h12 within each location pre vs post)
#meta$population <- with(meta, paste('Location')) #randomise by location only for h12

# Reorder and remove columns
column.order <- c('sample_id', 'Location', 'control_phase','population')
meta <- meta[, ..column.order]

# A function that will return a list of length k, where each k is a random shuffling of the input vector
shuffle <- function(x, k){
  data.table(replicate(k, sample(x, length(x), replace = F)))
}
# Create 1000 or more randomisations of the phenotype labels, stratified by population, and add them to the control
# table
set.seed(42)
num.randomisations <- 10000
replicate.names <- paste('r', str_pad(1:num.randomisations, nchar(as.character(num.randomisations)), pad = 0), sep = '')

meta[, (replicate.names) := shuffle(control_phase, num.randomisations), by = population] #use location for h12

# Write the table to file
fwrite(meta, './fst/data/location_randomisations_pop.csv', sep = '\t', quote = F)

