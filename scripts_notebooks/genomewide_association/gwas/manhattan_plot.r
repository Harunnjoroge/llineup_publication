############################################################
## Manhattan Plot with CRE Annotations 
############################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggfx)
library(ggrepel)

setwd('/home/harunnn/lstm_scratch/network_scratch/llineup/llineup_publication/Data/cre_enrichment/results')

############################################################
## 1. Load input data
############################################################

cre_map <- fread("SNP_CRE_mapping.tsv")
snp <- fread("/home/harunnn/lstm_scratch/network_scratch/llineup/llineup-genomics/glm/allnets_results/gwas_snps.csv")

############################################################
## 2. Merge CRE annotations
############################################################

merged <- snp %>%
  left_join(
    cre_map %>% select(CHR, BP, CRE_annotation),
    by = c("CHR", "BP")
  )

merged$CRE_annotation[is.na(merged$CRE_annotation)] <- "non_THS"

############################################################
## 3. Compute score metrics
############################################################

merged$abs_score    <- -log10(merged$p)
merged$signed_score <- sign(merged$coef) * merged$abs_score

############################################################
## 4. Force chromosome order
############################################################

chrom_order <- c("2R","2L","3R","3L","X")
merged$CHR <- factor(merged$CHR, levels = chrom_order, ordered = TRUE)

############################################################
## 5. Compute cumulative positions
############################################################

merged <- merged %>%
    arrange(CHR, BP) %>%
    group_by(CHR) %>%
    mutate(chr_len = max(BP)) %>%
    ungroup()

chr_sizes <- merged %>%
    group_by(CHR) %>%
    summarize(chr_len = max(BP))

chr_offsets <- chr_sizes %>%
    mutate(offset = lag(cumsum(chr_len), default = 0))

merged <- merged %>%
    left_join(chr_offsets, by = "CHR") %>%
    mutate(BP_cum = BP + offset)

############################################################
## 6. Top 20 SNPs
############################################################

top20 <- merged %>%
    arrange(desc(abs_score)) %>% 
    slice(1:20)

merged$is_top20 <- paste(merged$CHR, merged$BP) %in%
                   paste(top20$CHR, top20$BP)

############################################################
## 7. CRE colors
############################################################

cre_cols <- c(
    "Promoter"            = "#E41A1C",
    "5' UTR"              = "#389ff3ff",
    "3' UTR"              = "#def112ff",
    "Intronic"            = "#984EA3",
    "Distal Intergenic"   = "#f6be06ff",
    "Exonic"              = "#e35e12ff",
    "Downstream (<1kb)"   = "#F781BF",
    "Downstream (1-2kb)"  = "#1ebd2cff",
    "Downstream (2-3kb)"  = "#1e4d3eff",
    "non_THS"             = "#0c17eaff"
)

############################################################
## 8. Chromosome bands — fully fixed (no chr_len clashes)
############################################################

# Rename before joining to avoid suffixes
chr_sizes2 <- chr_sizes %>% rename(chr_len_size = chr_len)
chr_offsets2 <- chr_offsets %>% rename(chr_len_offset = chr_len)

bands <- chr_sizes2 %>%
    left_join(chr_offsets2, by="CHR") %>%
    mutate(
        xmin = offset,
        xmax = offset + chr_len_size,
        fill_id = row_number() %% 2,
        band_col = ifelse(fill_id == 0, "#F0F0F0", "#c0bdbdff")
    )

############################################################
## 9. Manhattan function using ggfx raster
############################################################

plot_manhattan <- function(df, yvar = "signed_score", title = "") {

    ylab <- ifelse(yvar == "signed_score",
                   "Signed -log10(p)",
                   "-log10(p)")

    df_bg  <- df %>% filter(!is_top20)
    df_top <- df %>% filter(is_top20)

    p <- ggplot() +

        # Chromosome background bands
        geom_rect(
            data = bands,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = bands$band_col,
            color = NA
        ) +

        # Rasterized background points
        with_raster(
            geom_point(
                data = df_bg,
                aes(x = BP_cum, y = .data[[yvar]]),
                size = 0.6,
                alpha = 0.7,
                color = "grey30"
            ),
            dpi = 900
        ) +

        # Vector top SNPs
        geom_point(
            data = df_top,
            aes(x = BP_cum, y = .data[[yvar]], fill = CRE_annotation),
            shape = 21,
            size = 2.2,
            color = "black",
            stroke = 0.22
        ) +

        scale_fill_manual(values = cre_cols) +

        scale_x_continuous(
            breaks = chr_offsets2$offset + chr_sizes2$chr_len_size/2,
            labels = chr_sizes2$CHR,
            expand = expansion(mult = c(0.01, 0.02))
        ) +
        labs(
            x = NULL,
            y = ylab,
            title = title
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        )

    return(p)
}

############################################################
## 10. Build and save plots
############################################################

p_signed <- plot_manhattan(
    merged, "signed_score", "Signed Manhattan Plot (Top 20 SNPs)"
)
p_unsigned <- plot_manhattan(
    merged, "abs_score", "Unsigned Manhattan Plot (Top 20 SNPs)"
)

ggsave("CRE_manhattan_signed_top20.svg",   p_signed,   width=12, height=5)
ggsave("CRE_manhattan_unsigned_top20.svg", p_unsigned, width=12, height=5)

message("✅ All Manhattan plots generated:")
message("  - CRE_manhattan_signed_top20.svg")
message("  - CRE_manhattan_unsigned_top20.svg")


~/lstm_scratch/network_scratch/gaard_kenya/gwas/scripts/pre_commit.sh