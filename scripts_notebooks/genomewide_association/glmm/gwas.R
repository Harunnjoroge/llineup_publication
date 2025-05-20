# Fit a mixed model to look at the association between intervention phase vs change in snps frequencies
#######################################################################################################################################
library(data.table)
library(ggplot2)
library(stringr)
library(glmmTMB)
library(grid)
library(gridExtra)
library(gtable)
library(tidyverse)
library(future.apply)

#working directory
wd <- "/llineup_publication/Data"
setwd(wd)

# Load the genotyping data
# fetch genotype subset  files 

dir.genotype_subset <- './genotype_subset' #create
genotype.files <- list.files(dir.genotype_subset)
genotype.files <- genotype.files [grepl('.*rds', genotype.files)]
#gwas round 1 vs 3
dir.glm.result <- './allnets_results'  #set the rounds you want 
#dir.glm.result <- './glm/glm_test'
#read subset genotype files  and run gwas
for (c in genotype.files){
  file.name <- sub('.rds','',c)
  setwd(dir.genotype_subset)
  gt  <- readRDS(c)
  setwd(wd)
  metadata <- fread('./merged_llineup_metadata.csv',na.strings = '') # load your metadata
  # rounds ilter 1 and 5
  sample.names <- colnames(gt)[3:ncol(gt)]
  metadata_allnets<- metadata[sample_id %in% sample.names & sex_call=='F' & RND %in% c(1,5)]
  
  #Prep data-Select pbo or npbo samples in round 1 and 3 or any net and rounds you want
 # metadata_pbo<- metadata[sample_id %in% sample.names & sex_call=='F' & RND %in% c(1,3)&LLIN.actual=='PBO LLIN']

  #metadata_npbo<- metadata[sample_id %in% sample.names & sex_call=='F' & RND %in% c(1,3)&LLIN.actual=='Non-PBO LLIN']
  metadata <- metadata_allnets %>% distinct(sample_id, .keep_all = TRUE)#Remove duplicated row in meta
  gt.id <- gt$`Unnamed: 0`
  #gt <- gt[,1:2 := NULL]
  gt <- gt %>% select(metadata$sample_id)
  setcolorder(gt, metadata$sample_id)
  
  # We will turn the gt to haplotypes
  gen2hap <- function(n){
    outvec <- rep(1,2)
    outvec[seq_len(2-n)] <- 0
    outvec
  }
  
  gen2hap.vector <- function(N){
    unlist(lapply(N, gen2hap))
  }
  
  #run in parallel
  # multi cores processing
  n.cores <- 35
  # Set the number of cores that future_apply will use for all its operations
  plan(tweak(multisession, workers = n.cores))
  
  #allow variables larger than default 500mb
  options(future.globals.maxSize=+Inf)
  #genotype to haplotype
  start.time <- Sys.time()
   cat('Haplotype conversion', c, 'started',start.time, '\n')
  ug.hap <- future_apply(gt, 1, gen2hap.vector) %>% data.frame()
  
  end.time <- Sys.time()

  cat('Haplotype conversion', c, 'done at',end.time, '\n')
  
  colnames(ug.hap) <- gt.id
  
  #double metadata because Turning genotype to haps duplicated each snp position
  
  model.meta<- metadata[,.(sample_id= rep(sample_id, each =2),
                           RND = rep(RND, each = 2),
                           HSD = rep(HSD, each = 2),
                           EA = rep(EA, each = 2),
                           HHID.rnd = rep(HHID.rnd, each = 2),
                           Wave = rep(Wave, each = 2),
                           LLIN.actual = rep(LLIN.actual, each = 2),
                           Location = rep(Location, each = 2),
                           wgs.sample.id = rep(wgs.sample.id, each = 2)
  )]
  
  #Running the model
  rnd<- as.factor(model.meta$RND)
  hsd<- as.factor(model.meta$HSD)
  ea<- as.factor(model.meta$EA)
  hhid <- as.factor(model.meta$HHID.rnd)
  wave <- as.factor(model.meta$Wave)
  location<- as.factor(model.meta$Location)
  llin <- as.factor(model.meta$LLIN.actual)
  
  #Extract p values in coefficients
  get.term.details <- function(model, test, model.term, test.term){
    p <- test[['Pr(>Chi)']][which(attributes(test)$row.names == test.term)]
    coefficient <- fixef(model)$cond[model.term]
    c(P = signif(p, 2), coef = signif(coefficient, 2))
  }
  
  #glmm model function
  run.model <- function(markers, ug.hap){
    
    model <- glmmTMB(get(markers) ~ rnd + location + llin + (1|hsd), 
                     data = ug.hap, family = 'binomial')
    test1 <- drop1(model, test = 'Chisq')
    #print(test1)
    
    test.summary = get.term.details(model,test1,'rnd5','rnd')
    test.summary
  }
  
  markers <- gt.id
  
  # multi cores processing
  n.cores <- 45
  plan(tweak(multisession, workers = n.cores))
  #allow variables larger than default 500mb
  options(future.globals.maxSize=+Inf)
  
  #Run model
  start.time <- Sys.time()
  model.test <- do.call(rbind,future_lapply(markers, run.model,ug.hap))
  model.test <- as.data.frame(model.test) 
  model.test$snp.id <- gt.id
  saveRDS(model.test,paste(dir.glm.result,'/',file.name,'_glm.rds',sep = ''))
  cat('Subset glm', c, 'saved as ',file.name, '_glm.rds','\n')
  end.time <- Sys.time()
  cat(as.character('Glm of', c, 'started at', start.time, 'and finished after',end.time, '\n'))
  

}

# converted_time <- as.POSIXct(start.time, origin = "1970-01-01", tz = "BST")
# converted_time