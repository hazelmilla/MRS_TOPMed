#######################################################
############### JHS EPIC v1 - data prep ###############
#######################################################

library(dplyr)
library(tidyverse)

inputdir = "/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/Input_files"
outputdir = "/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/Output_files"
setwd(inputdir)

# Prepare data for analysis
pheno <- read.csv("pd_pheno_data.csv", header = TRUE)
dim(pheno) # [1] 1691   21
# has sample_name & cell types
subjid <- read.csv("pheno_020524.csv", header = TRUE)
dim(subjid) # [1] 1691   36
# has subject ID, Sample_Name
mrs <- read.csv("MRS_JHS_EPICV1.csv", header = TRUE)
dim(mrs) # [1] 1709    4
colnames(mrs)[1] <- "Sample_Name"
# has sample_name
analysis <- read.csv("analysis1_2020_update_dbgap_id.csv", header = TRUE)
dim(analysis)# [1] 3406  206
# has subject_ID
marital <- read.csv("marriage_status_dbgap_id.csv", header = TRUE)
dim(marital)# [1] 3400  4
# has subject_ID
pcs <- read.csv("PC_data_by_location.csv", header = TRUE)
dim(pcs)# [1] 1728   38
# has SUBJID (JID), SUBJECT_ID

# Merge pheno & subjid to find correspondence between 
# sample name and Subject ID
dat <- merge(pheno, subjid, by = "Sample_Name")
# Find correspondence between subject ID & mrs 
dat <- merge(dat, mrs, by = "Sample_Name")
# merge PCs to mrs data using subject ID
dat <- merge(dat, PCs, by = "SUBJECT_ID")
dim(dat)
#[1] 1691   96

data <- dat %>% subset(select = c("SUBJECT_ID", "SUBJID", "Sample_Name", "PC1", 
                                  "PC2", "PC3", "PC4", "PC5", "PC6", "PC7",
                                  "PC8", "PC9", "PC10", "NK", "Mono", "Gran", 
                                  "Bcell", "CD8T", "CD4T", "MRS", "MRS_FDR",
                                  "MRS_bonferroni"))
dim(data)
# [1] 1691   22

med <- median(data$MRS)
data$MRS2cat[data$MRS <= med] <- 1
data$MRS2cat[data$MRS > med] <- 2
table(data$MRS2cat)
#1   2 
#846 845 

med_fdr <- median(data$MRS_FDR)
data$MRS2cat_FDR[data$MRS_FDR <= med_fdr] <- 1
data$MRS2cat_FDR[data$MRS_FDR > med_fdr] <- 2
table(data$MRS2cat_FDR)
#1   2 
#846 845 

med_bonferroni <- median(data$MRS_bonferroni)
data$MRS2cat_bonferroni[data$MRS_bonferroni <= med_bonferroni] <- 1
data$MRS2cat_bonferroni[data$MRS_bonferroni > med_bonferroni] <- 2
table(data$MRS2cat_bonferroni)
#1   2 
#846 845 

# Add stress & covariates using analysis data
glimpse(analysis)
df <- merge(analysis, marital, by = "SUBJECT_ID")

# Model: control for age, sex, BMI, smoking, income, education, marital status
# alcohol, neighborhood SES
# age, BMI - continuous
# sex, smoking, alcohol - categorical
# Income - 4 categories
# Education - 3 categories
# marital status - 2 categories
# nbSESanascore (neighborhood SES) - continuous z score

# individual SES
# Income
# Create graph with missing info
sum(is.na(df$sex)) # no missing data
sum(is.na(df$age)) # no missing data

#BMI
BMI_missing <- sum(is.na(df$BMI)==TRUE) #7
BMI_count <- sum(is.na(df$BMI)==FALSE) #3393

#Smoking
smoking_missing <- sum(is.na(df$currentSmoker)==TRUE) #25
smoking_count <- sum(is.na(df$currentSmoker)==FALSE) #3375

#income
income_missing <- sum(is.na(df$Income)==TRUE) #509
income_count <- sum(is.na(df$Income)==FALSE) #2891

# education
education_missing <- sum(is.na(df$edu3cat==TRUE)) #11
education_count <- sum(is.na(df$edu3cat)==FALSE) #3389

# marriage
marital_missing <- sum(is.na(df$SOCA1A==TRUE)) #0
marital_count <- sum(is.na(df$SOCA1A)==FALSE) #3400

#alcohol
alc_missing <- sum(is.na(df$alc)==TRUE) #19
alc_count <- sum(is.na(df$alc)==FALSE) #3381

# Census Tract SES Score (Diez-Roux 1990)
SES_DiezRoux_missing <- sum(is.na(df$nbSESanascore)==TRUE) #3
SES_DiezRoux_count <- sum(is.na(df$nbSESanascore)==FALSE) #3397

dat_info <- data.frame(BMI = c(BMI_missing, BMI_count),
                       Smoking = c(smoking_missing, smoking_count),
                       Alcohol = c(alc_missing, alc_count),
                       Income = c(income_missing, income_count),
                       Education = c(education_missing,education_count),
                       Marital = c(marital_missing, marital_count),
                       SES = c(SES_DiezRoux_missing, SES_DiezRoux_count),
                       row.names = c("Missing", "Count"))

dat_info <- dat_info %>% t() %>% as.data.frame() %>% 
  mutate(Total = Missing + Count) %>% t() %>% as.data.frame()

library(pixiedust)
library(kableExtra)

dust(dat_info, keep_rownames = TRUE) %>%
  sprinkle_colnames(.rownames = " ") %>%
  kable() %>% kable_styling()

variables <- subset(df, select = c("SUBJECT_ID", "sex", "age", "BMI", 
                                   "currentSmoker", "Income", "edu3cat", "SOCA1A", 
                                   "nbSESanascore", "alc", "perceivedStress"))

stress_median <- summary(variables$perceivedStress)[3]
variables$stress2cat[variables$perceivedStress <= stress_median] <- 1
variables$stress2cat[variables$perceivedStress > stress_median] <- 2
table(variables$stress2cat)
#1    2 
#1737 1637 

# More even split relative to 1 is < median/2 is >= median
# Distributions when calculations greater than/equal 
# 1    2 
# 1470 1904 

# Split into married/widowed/sep vs never married
variables$marital2cat[variables$SOCA1A == "N"] <- "Not_married"
variables$marital2cat[variables$SOCA1A == "M"] <- "Married"
variables$marital2cat[variables$SOCA1A == "D"] <- "Not_married"
variables$marital2cat[variables$SOCA1A == "S"] <- "Not_married"
variables$marital2cat[variables$SOCA1A == "W"] <- "Not_married"
table(variables$marital2cat)
#Married Never_married 
#2983           410 

dim(variables)
# [1] 3400   13
write.table(variables, "variables_covariates.csv", sep = ",", row.names = FALSE,
            col.names = TRUE)

data <- merge(data, variables, by = "SUBJECT_ID")
write.table(data, "prepared_data.csv", sep = ",", row.names = FALSE,
            col.names = TRUE)

###### Prepare CHD Data for Cox Regression ######
chd_dat <- read.csv("incevtchd_2020_update_dbgap_id.csv", 
                    header = TRUE)
analysis <- read.csv("analysis1_2020_update_dbgap_id.csv", 
                     header = TRUE)
mrs <- read.csv("prepared_data.csv", 
                     header = TRUE)

CHDHx <- filter(analysis, CHDHx == 1) # no of participants = 256
noCHDHx <- filter(analysis, CHDHx != 1) # no of participants = 3150
# Remove participants with prior history of CHD

# all visits are visit 1 in this dataset
noCHDHx <- mutate(noCHDHx, date = as.Date(noCHDHx$VisitDate, 
                                          format =  "%m/%d/%Y"))
# Convert the date columns to Date type
chd_dat <- chd_dat %>%
  mutate(
    date_format = as.Date(date, format = "%m/%d/%Y"),
    v1date_format = as.Date(V1date, format = "%m/%d/%Y")
  )

# Calculate the time difference in days
chd_dat1 <- chd_dat %>%
  mutate(
    time_to_event = difftime(date_format, v1date_format, units = "days")
  ) %>%
  mutate(
    days_to_event = as.numeric(time_to_event))

data <- merge(noCHDHx, chd_dat1, by = "SUBJECT_ID")
data_chd <- subset(data, select = c("SUBJECT_ID", "CHD", "days_to_event"))

dat <- merge(data_chd, mrs, by = "SUBJECT_ID")


write.csv(dat, "chd_data.csv")

#######################################################
############### JHS EPIC v2 - data prep ###############
#######################################################

##### Data prep for MRS analysis in EPICv2 Data #####

#########################################################
############# PREPARE DNA METHYLATION DATA ##############
#########################################################

###### EPICV2 Data ######
library(dplyr)
library(tidyverse)

v2_path = "/proj/azannas/projects/JHS_Batch2"
path = "/work/users/h/e/hemilla/FVIII_analysis"
path1 = "/work/users/h/e/hemilla/"

setwd(v2_path)
load("beta.noXY.no.impute.RData")

setwd(path1)
subjid <- read.csv("subjid_TOE.csv", header = TRUE)
analysis <- read.csv("analysis1_2020_update_dbgap_id.csv", header = TRUE)

setwd("/work/users/h/e/hemilla/FVIII_analysis")
f8 <- read.delim("FVIII_JHS_dbGaP.txt", header = TRUE, sep="\t")

dat <- merge(analysis, f8, by = "SUBJECT_ID")

colnames(beta.noXY) # TOE ID
rownames(beta.noXY) # CPG site

names <- sapply(strsplit(rownames(beta.noXY), "_"), `[`, 1)
rownames(beta.noXY) <- names
dim(beta.noXY) #[1] 860960   1689

# find correspondence between TOEID and SUBJECT_ID
beta <- beta.noXY %>% t() %>% as.data.frame()
beta$TOEID <- rownames(beta)

tail(colnames(beta))
colnames(subjid)

beta1 <- merge(beta, subjid, by = "TOEID")

dim(beta1) #[1] 1687 860962
rownames(beta1)
colnames(beta1)

rownames(beta1) <- beta1$SUBJECT_ID
rownames(beta1)
colnames(beta1)

drop <- c("TOEID", "SUBJECT_ID")
beta1 = beta1[,!(colnames(beta1) %in% drop)]

beta <- beta1 %>% as.matrix() %>% t()
dim(beta) #[1] 860960   1687
colnames(beta) # subjid
rownames(beta) # cpg

setwd(path)
save(beta, file = "beta.v2.subjid.RData")

###############################################
############# PREPARE PHENO DATA ##############
###############################################

library(dplyr)
library(tidyverse)

path = "/work/users/h/e/hemilla/WHI_JHS_MRS/EPICv2"
setwd(path)

data <- read.csv("analysis1_2020_update_dbgap_id.csv", header = TRUE)
pc <- read.delim("JHS_pcair_pcs_output_TOEID.txt", header = TRUE, sep = "\t")
horvath <- read.csv("DNAmAgeCalcProject_14998_Results.csv", header = TRUE)
marital <- read.csv("marriage_status_dbgap_id.csv")
subjid <- read.csv("subjid_TOE.csv", header = TRUE)
mrs <- read.csv("MRS_JHS_EPICV2.csv", header = TRUE)

# get pc & horvath outputs together for genetic ancestry pcs 
# and Houseman cell type proportions
colnames(horvath)[1] = "TOEID"
pc_horvath <- merge(pc, horvath, by = "TOEID")
pc_ctp <- merge(pc_horvath, subjid, by = "TOEID")

# merge with data and marital to have all phenotypes in 1 place
pheno <- merge(data, marital, by = "SUBJECT_ID")
dat <- merge(pheno, pc_ctp, by = "SUBJECT_ID")

# subset the data 
dat <- dat %>% subset(select = c("SUBJECT_ID", "perceivedStress", "sex", "age",
                                 "VisitDate", "alc", "currentSmoker", "BMI", 
                                 "Income", "edu3cat", "nbSESanascore",
                                 "SOCA1A", "PC1", "PC2", "PC3", "PC4", "PC5", 
                                 "PC6", "PC7", "PC8", "PC9", "PC10", "NK", 
                                 "Bcell", "Mono", "Gran", "CD8T", "CD4T"))

# merge with MRS scores
data <- merge(dat, mrs, by = "SUBJECT_ID")

med <- median(data$perceivedStress, na.rm = TRUE)

data$stress2cat[data$perceivedStress <= med] <- 1
data$stress2cat[data$perceivedStress > med] <- 2

med_mrs <- median(data$MRS, na.rm = TRUE)

data$MRS2cat[data$MRS <= med_mrs] <- 1
data$MRS2cat[data$MRS > med_mrs] <- 2

med_mrs_FDR <- median(data$MRS_FDR, na.rm = TRUE)

data$MRS2cat_FDR[data$MRS_FDR <= med_mrs_FDR] <- 1
data$MRS2cat_FDR[data$MRS_FDR > med_mrs_FDR] <- 2

med_mrs_bonf <- median(data$MRS_bonferroni, na.rm = TRUE)

data$MRS2cat_bonf[data$MRS_bonferroni <= med_mrs_bonf] <- 1
data$MRS2cat_bonf[data$MRS_bonferroni > med_mrs_bonf] <- 2

data$marital2cat[data$SOCA1A == "N"] <- "Not_married"
data$marital2cat[data$SOCA1A == "M"] <- "Married"
data$marital2cat[data$SOCA1A == "D"] <- "Not_married"
data$marital2cat[data$SOCA1A == "S"] <- "Not_married"
data$marital2cat[data$SOCA1A == "W"] <- "Not_married"
table(data$marital2cat)

# save file
write.csv(data, "prepared_data_v2.csv", row.names = FALSE)

###### CHD data prep #######
path = "/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/V2"
setwd(path)

data_unfiltered <- read.csv("prepared_data_v2.csv", header = TRUE)
dim(data_unfiltered)
#[1] 1617   36

analysis <- read.csv("analysis1_2020_update_dbgap_id.csv", 
                     header = TRUE)

chd <- read.csv("incevtchd_2020_update_dbgap_id.csv", header = TRUE)
dim(chd)
#[1] 3406   14

# Data prep
chd_status <- analysis %>% subset(select = c("SUBJECT_ID", "CHDHx"))
chd_outcome <- chd %>% subset(select = c("SUBJECT_ID", "CHD", "V1date", "date"))

CHDHx <- filter(chd_status, CHDHx == 1) # no of participants = 256
noCHDHx <- filter(chd_status, CHDHx != 1) # no of participants = 3150

# Remove participants with prior history of CHD
data1 <- data_unfiltered %>% filter(data_unfiltered$SUBJECT_ID %in% noCHDHx$SUBJECT_ID)
dim(data1)
#[1] 1499   36

#convert date formats
chd_outcome <- chd_outcome %>%
  mutate(
    date_format = as.Date(date, format = "%m/%d/%Y"),
    v1date_format = as.Date(V1date, format = "%m/%d/%Y")
  )

# Calculate the time difference in days
chd_outcome <- chd_outcome %>%
  mutate(
    time_to_event = difftime(date_format, v1date_format, units = "days")
  ) %>%
  mutate(
    days_to_event = as.numeric(time_to_event))

data_unfiltered <- merge(data1, chd_outcome)
write.csv(data_unfiltered, "chd_data_v2.csv", row.names = FALSE)

#######################################################
################### MESA data prep ####################
#######################################################

library(tidyverse)
library(dplyr)
library(data.table)
library(xml2)
library(purrr)
library(dplyr)

# Get covariates
path = "/work/users/h/e/hemilla/WHI_MESA_MRS"
setwd(path)
cov <- read.csv("Input/MESA_methyl_visit_1_cov.csv", header = TRUE) 
# contains technical covariates, genetic ancestry PCS, & cell type proportions
# also includes age

mrs <- read.csv("Input/MRS_MESA.csv", header = TRUE)

setwd("/proj/azannas/projects/MESA/Phenotypes MESA/CSV/SHAReE1_ZannasPheno_20230626")
phenos <- fread("SHAReE1_ZannasPheno_20230623.csv", header = TRUE)
setwd("/proj/azannas/projects/MESA/Phenotypes MESA/CSV/SHAReE1_ZannasEV_20230626")
chd <- fread("SHAReE1_ZannasEV_20230623.csv", header = TRUE)

setwd("/proj/azannas/projects/MESA/Phenotypes dbGaP/")
ses <- fread("phs000209.v13.pht004321.v1.p3.c1.MESA_AncilMesaNeighborSESExam1.HMB.txt", 
             skip = 10, header="auto", sep="auto", dec=".", stringsAsFactors = FALSE)


# phenotypes needed: sex, age, currentSmoker, alc, BMI, Income, education, 
# marital status, SES, genetic ancestry PCs, technical covariates, CT proportions
# for Cox regression: CHD & days to event (CHD), omit CHDHx (CHD history)

# stress phenotypes
stress_pheno <- subset(phenos, select = c("sidno", "hprb1pt1", "hprb2pt1", 
                                          "hprb3pt1", "hprb1ot1", "hprb2ot1", 
                                          "hprb3ot1", "job1prb1", "job2prb1", 
                                          "job3prb1", "mon1prb1", "mon2prb1", 
                                          "mon3prb1", "rel1prb1", "rel2prb1",
                                          "rel3prb1"))
# hprb - health problems (self)
# hprbot - health problems (someone close)
# jobprb - job problems
# monprb - financial strain
# relprb - relationship problems
# 1 - ongoing, 2 - ongoing for >6mo, 3 - level of stress (1-3)

stress_pheno <- stress_pheno %>% mutate(
  hprb = case_when(
    hprb2pt1 == 1 ~ hprb2pt1 + hprb3pt1,
    is.na(hprb2pt1) ~ NA_real_,
    TRUE ~ 0)
  ) %>% mutate(
    hprb_ot = case_when(
      hprb2ot1 == 1 ~ hprb2ot1 + hprb3ot1,
      is.na(hprb2pt1) ~ NA_real_,
      TRUE ~ 0)
  ) %>% mutate(
  jobprb = case_when(
    job2prb1 == 1 ~ job2prb1 + job3prb1,
    is.na(job2prb1) ~ NA_real_,
    TRUE ~ 0)
  ) %>% mutate(
  monprb = case_when(
    mon2prb1 == 1 ~ mon2prb1 + mon3prb1,
    is.na(mon2prb1) ~ NA_real_,
    TRUE ~ 0 )
  ) %>% mutate(relprb = case_when(
    rel2prb1 == 1 ~ rel2prb1 + rel3prb1,
    is.na(rel2prb1) ~ NA_real_,
    TRUE ~ 0))

stress_pheno <- mutate(stress_pheno, chronic_stress = hprb + hprb_ot + jobprb + monprb + relprb)
hist(stress_pheno$chronic_stress)

stress_med <- median(stress_pheno$chronic_stress, na.rm = TRUE)
stress_pheno$stress2cat[stress_pheno$chronic_stress < stress_med] <- 1
stress_pheno$stress2cat[stress_pheno$chronic_stress >= stress_med] <- 2
table(stress_pheno$stress2cat)

other_pheno <- subset(phenos, select = c("sidno", "gender1", "bmi1c", "bmicat1c", 
                                         "marital1", "income1", "educ1", 
                                         "curalc1", "cursmk1", "race1c"))

#sum(other_pheno$gender1 != other_pheno$phxsex1, na.rm = TRUE)
#[1] 0
#sum(is.na(other_pheno$gender1))
#[1] 0
#sum(is.na(other_pheno$phxsex1))
#[1] 22

other_pheno <- other_pheno %>% mutate( 
  marital2cat = case_when(
  marital1 == 1 ~ 1,
  marital1 %in% c(2, 3, 4, 5) ~ 0, 
  marital1 == 6 ~ NA_real_,
  is.na(marital1) ~ NA_real_)
) %>% mutate(
  edu4cat = case_when(
    educ1 %in% c(0, 1, 2) ~ 0,
    educ1 == 3 ~ 1,
    educ1 == 4 ~ 2, 
    educ1 %in% c(5, 6, 7, 8) ~ 3,
    is.na(educ1) ~ NA_real_)
) %>% mutate(
  income4cat = case_when(
    income1 %in% 1:8 ~ 0,
    income1 %in% 9:10 ~ 1, 
    income1 == 11 ~ 2,
    income1 %in% 12:13 ~ 3,
    is.na(income1) ~ NA_real_)
)
  
# 0 = <HS
# 1 = HS/GED
# 2 = >HS/GED

other_pheno$bmi3cat <- ntile(other_pheno$bmi1c, 3)
hist(other_pheno$edu4cat)
hist(other_pheno$income4cat)

chd_pheno <- subset(chd, select = c("sidno", "chda", "chdatt", "mi", "mitt", 
                                    "prebase", "exall"))
# chda - chd all 
# chdatt - time to chd (days)
# mi - myocardial infarction
# mitt - time to MI (days)
# prebase - type of prebase event
# exall - exclusion due to prebase event (1 = yes)

# Prepare MRS data
# Add categories for MRS data
mrs_med <- median(mrs$MRS, na.rm = TRUE)
mrs$MRS2cat[mrs$MRS <= mrs_med] <- 1
mrs$MRS2cat[mrs$MRS > mrs_med] <- 2
table(mrs$MRS2cat)

mrs_med_bonf <- median(mrs$MRS_bonf, na.rm = TRUE)
mrs$MRS2cat_bonf[mrs$MRS_bonf <= mrs_med_bonf] <- 1
mrs$MRS2cat_bonf[mrs$MRS_bonf > mrs_med_bonf] <- 2
table(mrs$MRS2cat_bonf)

colnames(mrs)[1] <- "TOEID"

colnames(cov)[2:3] <- c("sidno", "TOEID")

cov_pheno <- merge(other_pheno, cov, by = "sidno")
# merge technical covariates, ancestry PCs, and CT proportions with covariates
dat <- merge(cov_pheno, stress_pheno, by = "sidno") # merge with stress data
dat <- merge(dat, chd_pheno, by = "sidno") # merge with chd data
dat <- merge(dat, mrs, by = "TOEID") # merge with MRS data
dat <- merge(dat, ses, by = "sidno", all.x = TRUE) # merge with SES data
# higher value of F1_PC2_1 corresponds to lower SES?
dim(dat) #[1] 870  88

colSums(is.na(dat))

setwd(path)
write.csv(dat, "Input/prepared_data_MESA.csv", row.names = FALSE)

# genetic ancenstry PCs: Lowercase "pc"+"1-10"
# capitalized PCs are LV PCs




