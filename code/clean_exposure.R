# Load required packages

library(magrittr)

# Set p-value threshold for filtering GWAS results

p <- 1e-5

# For each file...

for (i in c("bmiannotated","waist_circumannotated","waist_hip_ratioannotated")) {
  
  # Load data
  
  df <- data.table::fread(paste0("./",i,".txt.gz"), data.table=FALSE)
  
  # Restrict to SNPs that meet the p-value threshold
  
  df <- df[as.numeric(df$pval_fe)<p & !is.na(df$pval_fe),]
  
  # Combine information on allele frequency and sample size
  
  df <- df %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(af = min(af_uganda, af_DCC, af_DDS, af_AADM, na.rm = TRUE),
                  n = sum(no_uganda, no_DCC, no_DDS, no_AADM, na.rm = TRUE))
  
  # Seperate SNP ID into component parts
  
  df <- tidyr::separate(df, snpid, into = c("chr","pos","ea","oa"), sep = ":")
  
  # Restrict to required columns
  
  df <- df[,c("chr","pos","ea","oa","beta_fe", "se_fe", "pval_fe", "af", "n")]
  
  # Save restricted dataset
  
  data.table::fwrite(df,paste0("./restricted_",i,".csv"), row.names = FALSE)
  
}