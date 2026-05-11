## Residualized MRS & CHD in the MESA cohort

###### FDR-adjusted sites #######

setwd("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/002 Zannas_lab/MRS TOPMed/MRS_WHI_MESA/")
data_unfiltered <- read.csv("Input/prepared_data_MESA.csv", header = TRUE)

library(tidyverse)
library(dplyr)
library(car)
library(rstatix)
library(broom)
library(emmeans)
library(ggpubr)
library(survival)
library(survminer)

var <- c("MRS")
cov <- c("gender1", "AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
         "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
         "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
         "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
         "Bcell", "CD8T", "CD4T")
cov_s <- c("AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
         "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
         "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
         "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
         "Bcell", "CD8T", "CD4T")

data_unfiltered <- data_unfiltered %>% filter(prebase != "mi")
dim(data_unfiltered) #[1] 869  94

data <- data_unfiltered[complete.cases(data_unfiltered[ , c(var, cov)]), ]
dim(data) #[1] 707  94

sum(data$chda == 1) #47

resid <- paste("MRS ~", paste(cov, collapse = " + "))
resid_formula <- as.formula(resid)

res_fit <- lm(resid_formula, data=data)

data$MRS_resid <- resid(res_fit)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
# 1 - 354, 2 - 353
mesa_fdr_all <- data # 5-52

cox_model <- coxph(Surv(chdatt, chda) ~ MRSresid2cat, data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

mesa_cox_fdr <- coefficients

write.csv(mesa_cox_fdr, file = "Output/MESA_cox_results_fdr.csv")

######## Bonferroni-adjusted sites ########

setwd("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/002 Zannas_lab/MRS TOPMed/MRS_WHI_MESA/")
data_unfiltered <- read.csv("Input/prepared_data_MESA.csv", header = TRUE)

library(tidyverse)
library(dplyr)
library(car)
library(rstatix)
library(broom)
library(emmeans)
library(ggpubr)
library(survival)
library(survminer)

var <- c("MRS_bonf")
cov <- c("gender1", "AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
         "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
         "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
         "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
         "Bcell", "CD8T", "CD4T")
cov_s <- c("AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
           "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
           "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
           "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
           "Bcell", "CD8T", "CD4T")

data_unfiltered <- data_unfiltered %>% filter(prebase != "mi")
dim(data_unfiltered) #[1] 869  94

data <- data_unfiltered[complete.cases(data_unfiltered[ , c(var, cov)]), ]
dim(data) #[1] 707  94

resid <- paste("MRS_bonf ~", paste(cov, collapse = " + "))
resid_formula <- as.formula(resid)

res_fit <- lm(resid_formula, data=data)

data$MRS_resid <- resid(res_fit)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
# 1 - 354, 2 - 353
mesa_bonf_all <- data

cox_model <- coxph(Surv(chdatt, chda) ~ MRSresid2cat, data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

mesa_cox_bonf <- coefficients

write.csv(mesa_cox_bonf, file = "Output/MESA_cox_results_bonf.csv")

######## TNF sites ########

library(tidyverse)
library(dplyr)
library(survival)
library(survminer)

setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_MESA")

data_unfiltered <- read.csv("Input/prepared_data_MESA.csv", header = TRUE)

setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_MESA")
tnf_mrs <- read.csv("Input/TNF_sites_MRS_MESA.csv", header = TRUE)
data_unfiltered <- read.csv("Input/prepared_data_MESA.csv", header = TRUE)

glimpse(tnf_mrs)
colnames(tnf_mrs)[1] <- "TOEID"
glimpse(data_unfiltered)
dim(data_unfiltered) #[1] 870  94

data_unfiltered <- merge(tnf_mrs, data_unfiltered, by = "TOEID")

data_unfiltered <- data_unfiltered %>% filter(prebase != "mi")
dim(data_unfiltered) #[1] 869  96

var <- c("MRS_TNF", "TNF_avg")
cov <- c("gender1", "AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
         "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
         "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
         "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
         "Bcell", "CD8T", "CD4T")

data <- data_unfiltered[complete.cases(data_unfiltered[ , c(var, cov)]), ]
dim(data) #[1] 707  96

sum(data$chda == 1) #47

resid <- paste("MRS_TNF ~", paste(cov, collapse = " + "))
resid_formula <- as.formula(resid)

res_fit <- lm(resid_formula, data=data)

data$MRS_resid <- resid(res_fit)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat) # 1 - 354,   2 - 353

model <- coxph(Surv(chdatt, chda) ~ MRSresid2cat, data = data)

test_ph <- cox.zph(model)
print(test_ph)

coefficients <- summary(model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

write.csv(coefficients, file = "Output/MESA_cox_results_tnf.csv")
