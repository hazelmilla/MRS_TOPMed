# Data prep
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
