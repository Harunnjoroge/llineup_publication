#Plot Manhattan plots to visualise snps

library(tidyverse)
library(qqman)
library(fdrtool)
#merge the pvalues
wd <- '~/lstm_scratch/network_scratch/llineup/llineup-genomics/glm'
setwd(wd)
# fetch npbo corr results files 
dir.glm.results <- './allnets_results'
glm.files <- list.files(dir.glm.results)
glm.files <- glm.files [grepl('.rds', glm.files)]
glm.list <- list()

setwd(dir.glm.results)

#merge all ug haplotype npbo files
for (g in glm.files){
  glm.list [[g]] <- readRDS(g)
}

df <- do.call(rbind, glm.list)

df <- na.omit(df)

p <- df$P
snp.id <-df$snp.id

#fdr

fdr <- data.frame(fdrtool(p, statistic = 'pvalue', plot = F))
#df <- data.frame(snp_id,p)
df.adjusted <- data.frame(snp.id,p=fdr$qval)
df <- data.frame(snp.id,p=fdr$pval)
#significant filter
sig_adj<- df.adjusted %>% filter(p<0.05)
sig<- df %>% filter(p<0.01)
sig <- arrange(sig_adj, p)
sig_adj <- arrange(sig_adj, p)
snpsOfInterest <- df_sig$snp.id



#Manhattan adjusted
df.adjusted[c('CHR', 'BP')] <- str_split_fixed(df.adjusted$snp.id, ':', 2)
df.adjusted <- df.adjusted %>% mutate(CHR=fct_recode(CHR, "1"="2R", "2"="2L", "3" ="3R","4"="3L", "5"="X"))
df_n <- as.data.frame(apply(df.adjusted[,-1], 2, as.numeric))
df_n$snps.id <- df.adjusted$snp.id
df_n<- df_n[,c(4,1,2,3)]
df_n <- na.omit(df_n)
png('all_nets_adjusted_manhattan_highlighted.png', height =7, width=9, units='in',res = 300)
manhattan(df_n, snp = 'snps.id',chr = 'CHR', bp='BP', p='p', chrlabs = c("2R","2L","3R","3L","X"),
          suggestiveline = F,genomewideline = F, highlight = snpsOfInterest)

dev.off()

#Manhattan unadjusted
df[c('CHR', 'BP')] <- str_split_fixed(df$snp.id, ':', 2)
df <- df %>% mutate(CHR=fct_recode(CHR, "1"="2R", "2"="2L", "3" ="3R","4"="3L", "5"="X"))
#df_n <- as.data.frame(apply(df[,-1], 2, as.numeric))
df_n <- as.data.frame(apply(df[,!names(df) %in% c("snp.id", "coef.rnd5", "P")], 2, as.numeric))
df_n$snps.id <- df$snp.id
df_n<- df_n[,c(4,1,2,3)]
df_n <- na.omit(df_n)
png('all_nets_unadjusted_20_manhattan_highlighted_all.png', height =7, width=9, units='in',res = 300)
manhattan(df_n, snp = 'snps.id',chr = 'CHR', bp='BP', p='fdr_p', chrlabs = c("2R","2L","3R","3L","X"),
          suggestiveline = F,genomewideline = F, highlight = snpsOfInterest)
dev.off()

