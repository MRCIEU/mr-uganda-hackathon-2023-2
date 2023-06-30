##group work.
###doing the analysis of two sample MR for teh Hackathon

##clearing the global data lists
rm(list=ls())
##setiing the work directories
setwd("C:/hanchard lab 2019/hanchard lab 2019/hanchard lab 2019/BCM 2019/Managing innovation for sustainable health course 2021/MRC Mendelian randomization/Teaching material/DAY 5/mr-uganda-hackathon-2023-2-main/")


###loading the libraries
library(plyr) 
library(TwoSampleMR)
library(devtools)
library(calibrate)
library(ggrepel)
library(ggthemes)
library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)
library(png)

##reading in the data
ao <- read.csv("raw/exposure.csv")
##checking the first ten lines of the output

head (ao)
##this should list the datacross check
##checking the size of the exposures
dim(ao)
View(ao)
##this gives 550 11 columns

##checking the SNPS taht meet the P value thresholds of signofcance, Genome wide 5Xe-8
length(which(ao$pval_fe<=5E-8))
#these are 2 snps

##doing the one that shows the entire data

##getting the signifcant 
significant_exposure <- ao[ao$pval_fe<=5e-07,]
##dimensions of those with P value less than -07
dim(significant_exposure)
#this leaves 32 11 

##formatting the data 
names(ao)
#check out the two sample MR function of format_data
##creating a SNP list for merging
SNPS<-significant_exposure$snp
SNPS
##the format function renames the data into the new dataset for the Colnames in MR base
significant_format<-format_data(
  significant_exposure,
  type = "exposure",
  snps = NULL, ### PUT THE CREATED LIST HERE SO THAT THE RESTRICTED SNPS are replaced in the dtaformat created
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "snp",
  beta_col = "beta_fe",
  se_col = "se_fe",
  eaf_col = "af",
  effect_allele_col = "ea",
  other_allele_col = "oa",
  pval_col = "pval_fe",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "n",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  z_col = "z",
  info_col = "info",
  chr_col = "chr",
  pos_col = "pos",
  log_pval = FALSE
)


###
##creating a SNP list for merging/filtering
SNPS<-significant_format$SNP
SNPS 


###########################################Looking at the outcome data

##looking the outcome data
outcome <- read.csv("raw/outcome.csv")
dim(outcome)
#chcking teh outcome
colnames(outcome)
##RUNNING the format data command to ensure that the data has same SNpS
significant_outcome<-format_data(
  outcome,
  type = "outcome",
  snps = SNPS, ### PUT THE CREATED LIST HERE SO THAT THE RESTRICTED SNPS are replaced in the dtaformat created
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "snp",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "af",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "n",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  z_col = "z",
  info_col = "info",
  chr_col = "chr",
  pos_col = "pos",
  log_pval = FALSE
)
##the allele frequency is missing in the data for the outcome
colnames(significant_outcome)
#checking the dimensions of the data
dim(significant_outcome)
nrow(significant_outcome)


##clumping ensures the independence of the data. the data for the exposure must be formated.
#clumping the data. Use the LD clumping
#clumping window is about the compuational burden and speed.
 clumped_exposure <- clump_data(
  significant_format,
  clump_kb = 10000,
  clump_r2 = 0.01,#R squared to be changed to 0.01
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "AFR" #3use an African reference set.
)
 ##sometimes you may need to downlaod the data for local files 

 #############################3harmonize the data
 ##the clumped data with 12 rows is now kept
 clumped_exposure <- read.csv("clumped_data.csv")
 dat <- harmonise_data(clumped_exposure, significant_outcome, action = 2) ##checks for which is the the
View(dat) 
#two rows were lost

##check the unique comaprison
##looking at the data dimensions
dim(dat)
#this now has 20 rows and 32 cloumns 
##removing the palindromes
palindromic_at<-subset(dat,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(dat,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(dat,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(dat,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")
dim(palindromic_at)
dim(palindromic_ta)
dim(palindromic_gc)
dim(palindromic_cg)
#two snps have at

##mr.keep shows whether the steps 
##view

View(palindromic_at)

##code to check which methods to use
TwoSampleMR::mr_method_list()


############################running the MR code.
##MR-Base R package to estimate the effects using the IVW, MR-Egger, weighted median and weighted mode methods. Have a look at the mr_method_list() function
mr_results <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_results
##adding the confidence intervals
mr_results$or<- exp(mr_results$b)
mr_results$lci<-exp(mr_results$b-1.96*mr_results$se)
mr_results$hci<-exp(mr_results$b+1.96*mr_results$se)
#runing the b[1] means for the  first value in beta column.
##qnorm(0.975) instead of 1.96
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
mr_results$lci <- exp(mr_results$lci*mr_results$sd)
mr_results$hci <- exp(mr_results$hci*mr_results$sd)
results<-cbind.data.frame(mr_results$outcome,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$exposure,mr_results$lci,mr_results$hci,mr_results$or, mr_results$sd )

##saving the results data
write.csv(results,"finalresults.cv")


##or you can use the qnorm(0.975)
##viewing the data.
table(results$`mr_results$method`)

##checking teh Unique command
#3the data outputs the logodds. We need to have it as odds ratios
##dimensiosn of the final dataset
dim(results)
#this gives 6 columns and 18 rows.

#######################limiting the data to the ivw
resultsivw<-results[results$`mr_results$method` == 'Inverse variance weighted' ,]
##this works well
resultsivw
#these are only six.
###making the confidence intervals and odds ratios

##the data needs to have the same usnits for ease of interpretation.
#3the SD of the trait should be the one we used. 
# Perform analysis

#res <- TwoSampleMR::mr(dat)

# Calculate confidence interval

#res$b_lci <- res$b - qnorm(0.975)*res$se
#res$b_uci <- res$b + qnorm(0.975)*res$se


######running the sentivity analyses
#Is there evidence of heterogeneity in the genetic effects?
het <- mr_heterogeneity(dat)
het
# the q VALUES measure the heterogenity. Above 3
##based on the P values, there is liimited evidence of heterogeneity.

#pleiotropy
pleio <- mr_pleiotropy_test(dat)
pleio
res_single <- mr_singlesnp(dat)
res_single

##if there is alot of NA, the MR Agger may not have run.. As close to zero is the gaol in regards to intercept.




