# Load exposure data

exposure <- data.table::fread("raw/exposure.csv", data.table = FALSE)

# Load outcome data

outcome <- data.table::fread("raw/outcome.csv", data.table = FALSE)

# Explore thresholds for instrument selection

table(exposure[exposure$pval_fe<5e-8,]$phenotype)
table(exposure[exposure$pval_fe<1e-7,]$phenotype)
table(exposure[exposure$pval_fe<5e-7,]$phenotype)

# Load exposure data

exposure <- read.csv("raw/exposure.csv")

# Set p-value threshold for instrument selection

p <- 5e-7

# Restrict exposure data

exposure <- exposure[exposure$pval_fe<p,]

# Format exposure data

exposure <- TwoSampleMR::format_data(exposure,
                                     type = "exposure",
                                     phenotype_col = "phenotype",
                                     snp_col = "snp",
                                     beta_col = "beta_fe",
                                     se_col = "se_fe",
                                     eaf_col = "af",
                                     effect_allele_col = "ea",
                                     other_allele_col = "oa",
                                     pval_col = "pval_fe",
                                     samplesize_col = "n",
                                     chr_col = "chr",
                                     pos_col = "pos")

# Format outcome data

outcome <- TwoSampleMR::format_data(outcome,
                                    type = "outcome",
                                    snps = exposure$SNP,
                                    phenotype_col = "phenotype",
                                    snp_col = "snp",
                                    beta_col = "beta",
                                    se_col = "se",
                                    effect_allele_col = "effect_allele",
                                    other_allele_col = "other_allele",
                                    pval_col = "p_value",
                                    samplesize_col = "n",
                                    chr_col = "chr",
                                    pos_col = "pos")

# Clump data

clumped <- TwoSampleMR::clump_data(exposure,
                                   clump_kb = 10000,
                                   clump_r2 = 0.01,
                                   clump_p1 = 1,
                                   clump_p2 = 1,
                                   pop = "AFR")

data.table::fwrite(clumped, "data/clumped_data.csv", row.names = FALSE)

# Harmonise data

dat <- TwoSampleMR::harmonise_data(exposure_dat = clumped,
                                   outcome_dat = outcome,
                                   action = 2)

# Perform analysis

res <- TwoSampleMR::mr(dat)

# Calculate confidence interval

res$b_lci <- res$b - qnorm(0.975)*res$se
res$b_uci <- res$b + qnorm(0.975)*res$se

# Convert to standard deviation scale

res$unit <- ""
res$unit <- ifelse(res$exposure=="bmi","SD",res$unit)
res$unit <- ifelse(res$exposure=="waist_hip_ratio","SD",res$unit)
res$unit <- ifelse(res$exposure=="waist_circum","cm",res$unit)

res$sd <- 1
res$sd <- ifelse(res$exposure=="bmi",3.83,res$sd)
res$sd <- ifelse(res$exposure=="waist_hip_ratio",0.16,res$sd)

# Convert to odds ratio scale

res$or <- exp(res$b*res$sd)
res$lci <- exp(res$b_lci*res$sd)
res$uci <- exp(res$b_uci*res$sd)

# Save results

data.table::fwrite(res, "output/results.csv", row.names = FALSE)