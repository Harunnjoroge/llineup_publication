library(data.table)
library(glmmTMB)
library(future.apply)
library(progress)
library(ggplot2)
### -----------------------------
### Load data
### -----------------------------
wd <- "~/lstm_scratch/network_scratch/llineup"
setwd(wd)

gt  <- readRDS("./llineup-genomics/glm/genotype_subset/gt_1.rds") #load one subset to get samples used in gwas analysis
#only necessary to extract our actual sample size automatically but sample size included below so you can skip the filter steps

metadata <- fread("./llineup_publication/Data/merged_llineup_metadata.csv", na.strings = "")

sample.names <- colnames(gt)[3:ncol(gt)]
#get pbo cohort  only as it is the smallest
meta.sub <- metadata[sample_id %in% sample.names & sex_call == "F" & llin_actual == 'PBO LLIN']

### Keep rounds 1 and 5 only
meta.sub <- meta.sub[RND %in% c(1,5)]

### Sample IDs
r1.ids <- meta.sub[RND == 1]$sample_id  #pbo 185
r5.ids <- meta.sub[RND == 5]$sample_id  #pbo 53

### Diploid sample sizes
r1.sample.size <- length(r1.ids) * 2
r5.sample.size <- length(r5.ids) * 2

### Number of SNPs
n.snps <- 9898581

### Precompute covariate vectors once (outside parallel code) to avoid error from data.table when you do it inside future
loc_vec  <- c(meta.sub[RND == 1]$Location,    meta.sub[RND == 5]$Location)
hsd_vec  <- c(meta.sub[RND == 1]$HSD,         meta.sub[RND == 5]$HSD)

# Expand to diploid level
loc_rep <- rep(loc_vec, each = 2)
hsd_rep <- rep(hsd_vec, each = 2)

### -----------------------------
### GLMM functions
### -----------------------------
extract_pvalue <- function(model, test, model.term, test.term){
    p <- test[["Pr(>Chi)"]][which(rownames(test) == test.term)]
    return(p)
}
analysis.function<- function(r1.allele.count, r5.allele.count){
    r1.alleles <- c(rep(0, r1.allele.count), rep(1, r1.sample.size - r1.allele.count))
	r5.alleles <- c(rep(0, r5.allele.count), rep(1, r5.sample.size - r5.allele.count))
    

        # Use precomputed covariates; do NOT call data.table inside future functions
    df <- data.frame(
        genotype = c(r1.alleles, r5.alleles),
        RND      = c(rep("r1", r1.sample.size), rep("r5", r5.sample.size)),
        Location = loc_rep,
        HSD      = hsd_rep,
        stringsAsFactors = TRUE
    )

    model <- glmmTMB(genotype ~ RND + Location + (1|HSD),  # H3 test add LLIN
                     data=df, family="binomial")

    test <- drop1(model, test="Chisq")

    extract_pvalue(model, test,
                   model.term="RND5",
                   test.term="RND")
}


### -----------------------------
### Simulation function
### -----------------------------
simulation_H3 <- function(e=2, n.affected.snps=10, fdr.threshold=0.05){
  
    unaffected.p <- runif(n.snps - n.affected.snps)

    r1.freqs  <- runif(n.affected.snps)
    r1.counts <- round(r1.freqs * r1.sample.size)

    r5.freqs  <- r1.freqs / e
    r5.counts <- round(r5.freqs * r5.sample.size)

    affected.p <- mapply(
        FUN = analysis.function,
        r1.counts,
        r5.counts
    )

    pvals <- c(affected.p, unaffected.p)
    fdr   <- p.adjust(pvals, method="BH")

    true.positive.rate <- sum(fdr[1:n.affected.snps] < fdr.threshold)/n.affected.snps
    true.positive.rate 
}


set.seed(42)
n.cores <- 30
plan(tweak(multisession, workers = n.cores))
options(future.globals.maxSize=+Inf)



###------------------------------
### test several e for which gives 80% power
###------------------------------
run_simulation_H3 <- function(k, e, n.affected.snps){
    results <- future_sapply(
        1:k,
        function(i) simulation_H3(e = e, n.affected.snps = n.affected.snps),
        future.seed = TRUE
    )
    
    list(
        results = results,
        mean_true_positive = mean(results, na.rm = TRUE),
        ci = quantile(results, probs = c(0.025, 0.975), na.rm = TRUE)
    )
}

find_e_for_power_logging <- function(
    target_power = 0.8,
    e_start = 1,
    e_step = 0.1,
    e_max = 10,
    k = 500,
    n.affected.snps = 10,
    log_file = "e_power_log.csv",
    plot_file = "power_curve.svg"
){
    # Initialize log file if missing
    if (!file.exists(log_file)) {
        write_csv(
            tibble(
                e = numeric(),
                power = numeric(),
                ci_low = numeric(),
                ci_high = numeric()
            ),
            log_file
        )
    }
    
    e <- e_start
    
    repeat {
        message("\n====================================")
        message(" Testing e = ", e)
        message("====================================")
        
        # Run simulation 
        out <- run_simulation_H3(
            k = k,
            e = e,
            n.affected.snps = n.affected.snps
        )
        
        power   <- out$mean_true_positive
        ci_low  <- out$ci[1]
        ci_high <- out$ci[2]
        
        message(" → Mean power = ", round(power, 4))
        message(" → 95% CI = [", round(ci_low, 4), ", ", round(ci_high, 4), "]")
        
        # Append results safely 
        new_row <- tibble(
            e = e,
            power = power,
            ci_low = ci_low,
            ci_high = ci_high
        )
        write_csv(new_row, log_file, append = TRUE)
        
        # Load the log every loop
        history <- read_csv(log_file, show_col_types = FALSE)
        
        # Plot with CI ribbon
        p <- ggplot(history, aes(x = e, y = power)) +
            geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2) +
            geom_line(size = 1.2) +
            geom_point(size = 3) +
            geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
            labs(
                title = "Power vs Effect Size (e)",
                subtitle = "Mean Power with 95% CI",
                x = "Effect size (e)",
                y = "Power"
            ) +
            theme_minimal(base_size = 16)
        
        ggsave(plot_file, p, width = 8, height = 6)
        message(" → Updated plot saved to: ", plot_file)
        
        # Stopping condition
        if (power >= target_power) {
            message("\n*** Found e achieving ≥ ", target_power*100, "% power ***")
            message("*** e = ", e, " ***\n")
            return(list(e = e, result = out, log = history))
        }
        
        e <- e + e_step
        
        if (e > e_max) stop("Did not reach target power before hitting e_max.")
    }
}

res <- find_e_for_power_logging(
    target_power = 0.8,
    e_start = 1,
    e_step = 0.1,
    k = 500,
    n.affected.snps = 10,
    log_file = "./llineup_publication/scripts_notebooks/genomewide_association/simulations/e_power_log.csv",
    plot_file = "./llineup_publication/scripts_notebooks/genomewide_association/simulations/power_curve.svg"
)
