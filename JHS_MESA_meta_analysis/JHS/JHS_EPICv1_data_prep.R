#### JHS EPIC v1 - MRS data prep ######

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
