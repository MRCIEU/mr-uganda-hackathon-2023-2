# Create exposure dataset

exposure <- data.frame(chr = character(),
                       pos = character())

for (i in c("bmiannotated","waist_circumannotated","waist_hip_ratioannotated")) {
  
  tmp <- data.table::fread(paste0("raw/restricted_",i,".csv"), data.table=FALSE)
  
  tmp$phenotype <- gsub("annotated","",i)
  
  exposure <- rbind(exposure, tmp)
  
}  

exposure <- exposure[!is.na(exposure$beta_fe),]

# Create outcome dataset

outcome <- data.frame(chr = character(),
                       pos = character())

for (i in c("GCST90091238_buildGRCh37","GCST90091239_buildGRCh37")) {
  
  tmp <- data.table::fread(paste0("raw/restricted_",i,".csv"), data.table=FALSE)
  
  tmp$phenotype <- ""
  tmp$phenotype <- ifelse(i=="GCST90091238_buildGRCh37", "Smoking initiation (ever vs never)", tmp$phenotype)
  tmp$phenotype <- ifelse(i=="GCST90091239_buildGRCh37", "Smoking cessation (current vs former)", tmp$phenotype)
  
  tmp$n <- NA
  tmp$n <- ifelse(i=="GCST90091238_buildGRCh37", 10558, tmp$n)
  tmp$n<- ifelse(i=="GCST90091239_buildGRCh37", 4257, tmp$n)
  
  outcome <- rbind(outcome, tmp)
  
}  

outcome <- outcome[!is.na(outcome$odds_ratio),]

outcome <- dplyr::rename(outcome,
                         "chr" = "chromosome",
                         "pos" = "base_pair_location")

# Calculate standard errors for outcome data

outcome$beta <- log(outcome$odds_ratio)
outcome$beta_lci <- log(outcome$ci_lower)
outcome$beta_uci <- log(outcome$ci_upper)

outcome$se <- ((outcome$beta_uci-outcome$beta_lci)/qnorm(0.975))

outcome[,c("odds_ratio","ci_lower","ci_upper","beta_lci","beta_uci")] <- NULL

# Load reference file that contains chr, pos and rsID

ref <- data.table::fread("raw/AFR.bim",
                         data.table = FALSE)

ref <- ref[,c("V1","V4","V2")]
colnames(ref) <- c("chr","pos","snp")

# Annotate exposure and outcome with rsIDs

exposure <- merge(exposure, ref,
             by = c("chr","pos"))

outcome <- merge(outcome, ref,
                  by = c("chr","pos"))

# Save datasets

data.table::fwrite(exposure, "raw/exposure.csv", row.names = FALSE)
data.table::fwrite(outcome, "raw/outcome.csv", row.names = FALSE)