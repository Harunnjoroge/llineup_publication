# Associate snps with haplotype clusters to find those strongly associated with the clusters hence useful as tags for the swept haplotype
library(data.table)
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)

setwd('~/lstm_scratch/network_scratch/llineup/llineup-genomics')
# Load data
meta <- fread("./llineup_publication/Data/merged_llineup_metadata.csv")
geno <- fread("./data/haplotype_genotypes/gt_2R:28473444-28509726.csv")
geno_trimmed <- fread("./data/haplotype_genotypes/gt_trimmed_2R:28473444-28509726.csv")

# Transpose the dataframe to have SNPs as columns
transpose_df <- function(snp) {
  df_t <- data.frame(t(snp), stringsAsFactors = FALSE)
  colnames(df_t) <- df_t[1, ]
  df_t <- df_t[-1, , drop = FALSE]
  df_t <- cbind(sample_id = rownames(df_t), df_t)
  rownames(df_t) <- NULL
  return(df_t)
}

snps <- transpose_df(geno)
snps_haps <- transpose_df(geno_trimmed)

# Merge metadata with genotypic data
meta <- na.omit(merge(meta, snps, by = "sample_id", all = TRUE))

# Load haplotype cluster data lstm_scratch/network_scratch/llineup/llineup-genomics
haplotypes <- fread('./llineup_publication/scripts_notebooks/haplotype_clustering/cluster_assignments_2R:28463444-28499726.csv', key = 'sample_name')
#setnames(meta, 'sample_id', 'sample_name')

haplotypes <- haplotypes[, c("sample_id", names(haplotypes)[8:11]), with = FALSE]
haplotypes[, sample_name := as.character(sample_name)]
meta[, sample_name := sub("\\.0*$", "", as.character(sample_name))]
meta <- merge(meta, haplotypes, by = "sample_id", all.x = TRUE)

# Identify SNPs and clusters
snp_columns <- grep("^2R:", names(meta), value = TRUE)   # Adjust chromosome label as needed
clusters <- grep("^cluster_", names(meta), value = TRUE)

# Initialize a dataframe to store p-values
p_values_df <- data.frame(matrix(ncol = length(clusters), nrow = length(snp_columns)))
rownames(p_values_df) <- snp_columns
colnames(p_values_df) <- clusters

# Compute chi-squared tests for SNP-cluster associations
for (cluster in clusters) {
  for (snp in snp_columns) {
    contingency_table <- table(meta[[snp]], meta[[cluster]])
    chisq_result <- chisq.test(contingency_table)
    p_values_df[snp, cluster] <- chisq_result$p.value
  }
}

# Convert p-values to -log10 scale
p_values_df[p_values_df == 0] <- 1e-320
log_p_values_df <- -log10(p_values_df)

# Reshape for ggplot2
log_p_values_long <- melt(as.matrix(log_p_values_df))

# Extract numeric SNP positions
log_p_values_long <- log_p_values_long %>%
  mutate(SNP_position = as.numeric(gsub(".*:(\\d+)", "\\1", Var1))) %>%
  arrange(SNP_position) %>%
  mutate(Var1 = factor(Var1, levels = unique(Var1)))

# Select the top 10 most significant SNPs per cluster
top_5_per_cluster <- log_p_values_long %>%
  group_by(Var2) %>%
  slice_max(order_by = value, n = 10, with_ties = FALSE) %>%
  ungroup()

# Filter data for heatmap
top_5_snps <- top_5_per_cluster$Var1
log_p_values_long_5 <- log_p_values_long %>%
  filter(Var1 %in% top_5_snps)

log_p_values_long_plot_all_filtered <- log_p_values_long %>%
  filter(Var1 %in% log_p_values_long_5$Var1)


# Function to plot heatmap
plot_heatmap <- function(df, title) {
  ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "grey50",
                         limits = c(0, max(df$value, na.rm = TRUE))) +
    labs(x = "Clusters", y = "SNPs", fill = "-log10(p-value)") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    ggtitle(title)
}

# Save genotype-based heatmap
p1 <- plot_heatmap(log_p_values_long_plot_all_filtered , "Top 5 SNP-Cluster Associations (Genotypes)")
ggsave("./Notebooks/plots/dgk_heat_20.svg", plot = p1, width = 12, height = 10, dpi = 300)

# Filter top SNPs present in haplotypes
snps_in_haps <- grep("^X:", names(snps_haps), value = TRUE)
top_hap_snps_per_cluster <- top_5_per_cluster %>%
  filter(Var1 %in% snps_in_haps)

log_p_values_long_hap_filtered <- log_p_values_long %>%
  filter(Var1 %in% top_hap_snps_per_cluster$Var1)


# Create and save haplotype-based heatmap
if (nrow(log_p_values_long_hap_filtered) > 0) {
  p2 <- plot_heatmap(log_p_values_long_hap_filtered, "Top 5 SNP-Cluster Associations (Haplotypes)")
  ggsave("./Notebooks/plots/dgk_haps_heat_20.svg", plot = p2, width = 12, height = 10, dpi = 300)
}

# Function to find the top SNP per cluster
top_sig_snps <- function(df) {
  df %>%
    group_by(Var2) %>%
    slice_max(order_by = value, n = 20, with_ties = FALSE)

}

# Print top SNPs for genotypes
print("Top SNPs for Genotypes:")
df <- top_sig_snps(log_p_values_long)
fwrite(df, "./data/haplotype_genotypes/top_20_snps_genotypes_cyp9k1.csv")

# Print top SNPs for haplotypes
print("Top SNPs for Haplotypes:")
df_haps <- top_sig_snps(top_hap_snps_per_cluster)
fwrite(df_haps, "./data/haplotype_genotypes/top_snps_20_haps_cyp9k1.csv")

df_haps_c2 <- df_haps[df_haps$Var2 == 'cluster_2',]



##################################################################################################################################################
#Plot association between top genotype snp vs cluster

# Scatter plot with ggplot2 and fit line

plot_scatter <- function(data, snp_column, cluster_column) {
  ggplot(data, aes(x = factor(get(snp_column)), y = factor(get(cluster_column)))) +
    geom_point() +  # Scatter plot points
    #geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear regression line (no confidence interval)
    labs(title = paste("Correlation between SNP Genotype (", snp_column, ") and Haplotype Cluster", sep = ""),
         x = "SNP Genotype", 
         y = "Haplotype Cluster") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)) +
    theme(legend.position = "none") 
}

# Assuming your data frame is called `meta` and you want to plot '2R:28498747' vs 'cluster_1'
p <- plot_scatter(meta, '2R:28498747', 'cluster_1')

# Save the plot to a file
ggsave("./Notebooks/plots/trial_scatter.png", plot = p, width = 12, height = 10, dpi = 300)


