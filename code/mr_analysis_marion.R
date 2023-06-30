#Get working directory
getwd()

#Load libraries 
library(TwoSampleMR)
library(png)

#Load exposure and outcone data
exposure <- read.delim("./raw/exposure.csv", header = TRUE, sep = ",")

outcome <- read.delim("./raw/outcome.csv", header = TRUE, sep = ",")

dim(exposure)
dim(outcome)
colnames(exposure)
colnames(outcome)

#Subset valid SNPs
significant_exposure <- exposure[exposure$pval_fe<=5e-07,]
dim(significant_exposure)

exposure_data <- format_data(
  significant_exposure,
  type = "exposure",
  snps = NULL,
  header = TRUE,
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
  pos_col = "pos"
)

dim(exposure_data)

# Clumped exposure data (use small clump window since SNPS have dense LD)

cl_exposure_data <- clump_data(
  exposure_data,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "AFR"
)

cl_exposure_data <- read.delim("C:/Users/Marion/Downloads/clumped_data.csv", sep = ",", header = T)
dim(cl_exposure_data)

# Get snp IDs in exposure

snp_instruments <- cl_exposure_data$SNP

# Outcome data format
outcome_data <- format_data(
  outcome,
  type = "outcome",
  snps = snp_instruments,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "snp",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  samplesize_col = "n",
  chr_col = "chr",
  pos_col = "pos"
)

dim(outcome_data)

# Data clumping

# Harmonize 
dat <- harmonise_data(cl_exposure_data, outcome_data, action = 2)
colnames(dat)
dim(dat)

# Estimate the causal effects of exposures on outcomes
mr_results <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_results
results<-cbind.data.frame(mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval)
colnames(results)

# Calculate confidence interval

mr_results$b_lci <- mr_results$b - qnorm(0.975)*mr_results$se
mr_results$b_uci <- mr_results$b + qnorm(0.975)*mr_results$se

# Convert to standard deviation scale

mr_results$unit <- ""
mr_results$unit <- ifelse(mr_results$exposure=="bmi","SD",mr_results$unit)
mr_results$unit <- ifelse(mr_results$exposure=="waist_hip_ratio","SD",mr_results$unit)
mr_results$unit <- ifelse(mr_results$exposure=="waist_circum","cm",mr_results$unit)

mr_results$sd <- 1
mr_results$sd <- ifelse(mr_results$exposure=="bmi",3.83,mr_results$sd)
mr_results$sd <- ifelse(mr_results$exposure=="waist_hip_ratio",0.16,mr_results$sd)

# Convert to odds ratio scale

mr_results$or <- exp(mr_results$b*mr_results$sd)
mr_results$lci <- exp(mr_results$b_lci*mr_results$sd)
mr_results$uci <- exp(mr_results$b_uci*mr_results$sd)

#Export results
write.csv(mr_results,"mr_results.csv")

# Sensitivity analyses
# Heteroigeneity test
het <- mr_heterogeneity(dat)
colnames(het)

# Horizontal pleiotropy
pleio <- mr_pleiotropy_test(dat)
pleio
res_single <- mr_singlesnp(dat)
res_single

# Visualize causal effects of exposures on outcomes
png("scatter.png")
mr_scatter_plot(mr_results, dat)
dev.off()

# Generate a forest plot of each of the SNP effects
png("./forest.png")
mr_forest_plot(res_single)
dev.off()

