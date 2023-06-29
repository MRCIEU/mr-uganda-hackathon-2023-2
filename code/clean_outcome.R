# Identify list of SNPs to extract from outcome data

snplist <- data.frame(chr = character(),
                      pos = character())

for (i in c("bmiannotated","waist_circumannotated","waist_hip_ratioannotated")) {

  tmp <- data.table::fread(paste0("./restricted_",i,".csv"), data.table=FALSE)
  
  tmp <- tmp[,c("chr","pos")]
  
  snplist <- rbind(snplist, tmp)
  
}  

snplist <- unique(snplist)

# For each file...

for (i in c("GCST90091238_buildGRCh37","GCST90091239_buildGRCh37")) {
  
  # Load data
  
  df <- data.table::fread(paste0("./",i,".tsv.gz"), data.table=FALSE)
  
  # Restrict to SNPs that are in the restricted exposure data files
  
  df <- merge(df, 
              snplist, 
              all.y = TRUE, 
              by.x = c("chromosome","base_pair_location"),
              by.y = c("chr","pos"))
  
  # Restrict to required columns
  
  df <- df[,c("chromosome","base_pair_location",
              "effect_allele","other_allele",
              "odds_ratio","ci_lower","ci_upper",
              "p_value")]
  
  # Save restricted dataset
  
  data.table::fwrite(df,paste0("./restricted_",i,".csv"), row.names = FALSE)

}
