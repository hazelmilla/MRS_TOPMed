########################################
###### JHS EPICv1 Cox regression #######
########################################

# Get CHD data and load packages
dir = "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/002 Zannas_lab/MRS TOPMed/MRS_WHI_JHS/V1/"
setwd(dir)

data_unfiltered <- read.csv("Input_files/chd_data.csv", header = TRUE)
dim(data_unfiltered)
#[1] 1554   40

library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(broom)
library(emmeans)
library(survminer)
library(survival)

######## THRESHOLD 1 (FDR) ########
# Model - sex, age, smoking, bmi, income, education, alc, marital status,
# SES, PC1-10, Houseman cell types
variables_JHS <- c("MRS", "sex", "age", "currentSmoker", "BMI", "Income",
                   "edu3cat", "alc", "marital2cat", "nbSESanascore",
                   "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",
                   "PC9", "PC10", "NK", "Mono", "Gran", "Bcell", "CD8T",
                   "CD4T") # create list of variables to remove NAs

# Remove NAs
data <- data_unfiltered[complete.cases(data_unfiltered[ , variables_JHS]), ]
dim(data)
# [1] 1321   40
sum(data$CHD == 1, na.rm = TRUE) #97

# PLOT RESIDUALS
model_resid <- lm(data$MRS ~ sex + age + currentSmoker + BMI +
                    Income + edu3cat + alc + marital2cat + nbSESanascore +
                    PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + NK + Mono + Gran + 
                    Bcell + CD8T + CD4T, data=data)

data$MRS_resid <- resid(model_resid)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
#1 - 709    2 - 708
v1_fdr_all <- data

cox_model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat, 
                   data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")
v1_fdr_cox <- coefficients

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")
write.csv(coefficients, file = "v1_cox_results_fdr.csv")

######### THRESHOLD 2 (Bonferroni) ########

# Model - sex, age, smoking, bmi, income, education, alc, marital status,
# SES, PC1-10, Houseman cell types
variables_JHS <- c("MRS_bonferroni", "sex", "age", "currentSmoker", "BMI", "Income",
                   "edu3cat", "alc", "marital2cat", "nbSESanascore",
                   "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",
                   "PC9", "PC10", "NK", "Mono", "Gran", "Bcell", "CD8T",
                   "CD4T") # create list of variables to remove NAs
# Remove NAs
data <- data_unfiltered[complete.cases(data_unfiltered[ , variables_JHS]), ]
dim(data)
# [1] 1321   40

# PLOT RESIDUALS
model_resid <- lm(data$MRS_bonferroni ~ sex + age + currentSmoker + BMI +
                    Income + edu3cat + alc + marital2cat + nbSESanascore +
                    PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + NK + Mono + Gran + 
                    Bcell + CD8T + CD4T, data=data)

data$MRS_resid <- resid(model_resid)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
#1 - 709    2 - 708
v1_bonf_all <- data

cox_model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat, 
                   data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")
v1_bonf_cox <- coefficients

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")
write.csv(v1_bonf_cox, file = "v1_cox_results_bonf.csv")

######## MRS_TNF #########
library(dplyr)
library(tidyverse)
library(survival)
library(survminer)

setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/WHI_TNF_sites_JHS")

chd <- read.csv("Input/chd_data.csv", header = TRUE)
mrs <- read.csv("Input/TNF_sites_MRS_JHSv1.csv", header = TRUE)

colnames(mrs)[1] <- "Sample_Name"
data_unfiltered <- merge(mrs, chd, by = "Sample_Name")
data_unfiltered <- subset(data_unfiltered, select = -X)
dim(data_unfiltered)
#[1] 1554   41

data <- data_unfiltered %>% drop_na(CHD) %>% drop_na(days_to_event) %>%
  drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>% drop_na(TNF_MRS)
dim(data)
#[1] 1291   41

res_fit <- lm(TNF_MRS ~ sex + age + 
                currentSmoker + BMI + Income + edu3cat + alc + 
                marital2cat + nbSESanascore + 
                PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                PC8 + PC9 + PC10 + NK + Mono + Gran + 
                Bcell + CD8T + CD4T, data=data)

data$MRS_resid <- resid(res_fit)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat) # 1 - 646,   2 - 645

model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat, data = data)

test_ph <- cox.zph(model)
print(test_ph)

# Summarize the model
coefficients <- summary(model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

write.csv(coefficients, file = "Output/TNF_MRS_CHD_JHSv1.csv")

# Run Cox regression for average TNF methylation
res_fit <- lm(TNF_avg ~ sex + age + 
                currentSmoker + BMI + Income + edu3cat + alc + 
                marital2cat + nbSESanascore + 
                PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                PC8 + PC9 + PC10 + NK + Mono + Gran + 
                Bcell + CD8T + CD4T, data=data)

data$avg_resid <- resid(res_fit)

avg_median <- summary(data$avg_resid)[3]
data$avg.resid2cat[data$avg_resid <= avg_median] <- 1
data$avg.resid2cat[data$avg_resid > avg_median] <- 2
table(data$avg.resid2cat) # 1 - 646,   2 - 645

v1_data <- data

model <- coxph(Surv(days_to_event, CHD) ~ avg.resid2cat, data = data)

test_ph <- cox.zph(model)
print(test_ph)

# Summarize the model
coefficients <- summary(model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

write.csv(coefficients, file = "Output/TNF_avg_CHD_JHSv1.csv")

########################################
###### JHS EPICv2 Cox regression #######
########################################

dir = "/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/V2"
setwd(dir)

data_unfiltered <- read.csv("chd_data_v2.csv", header = TRUE)
dim(data_unfiltered)
#[1] 1499   43

library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(broom)
library(emmeans)
library(survminer)
library(survival)

######## THRESHOLD 1 (FDR) ########
# Model - sex, age, smoking, bmi, income, education, alc, marital status,
# SES, PC1-10, Houseman cell types
variables_JHS <- c("MRS", "sex", "age", "currentSmoker", "BMI", "Income",
                   "edu3cat", "alc", "marital2cat", "nbSESanascore",
                   "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",
                   "PC9", "PC10", "NK", "Mono", "Gran", "Bcell", "CD8T",
                   "CD4T") # create list of variables to remove NAs
data <- data_unfiltered[complete.cases(data_unfiltered[ , variables_JHS]), ]
dim(data)
# [1] 1227   43
sum(data$CHD == 1, na.rm = TRUE) #70

# PLOT RESIDUALS
model_resid <- lm(data$MRS ~ sex + age + currentSmoker + BMI +
                    Income + edu3cat + alc + marital2cat + nbSESanascore +
                    PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + NK + Mono + Gran + 
                    Bcell + CD8T + CD4T, data=data)

data$MRS_resid <- resid(model_resid)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
#1 - 614    2 - 613
v2_fdr_all <- data

cox_model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat, 
                   data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

v2_fdr_cox <- coefficients

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")
write.csv(v2_fdr_cox, file = "v2_cox_results_fdr.csv")


######## THRESHOLD 2 (BONFERRONI) ########
# Model - sex, age, smoking, bmi, income, education, alc, marital status,
# SES, PC1-10, Houseman cell types
variables_JHS <- c("MRS_bonferroni", "sex", "age", "currentSmoker", "BMI", "Income",
                   "edu3cat", "alc", "marital2cat", "nbSESanascore",
                   "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",
                   "PC9", "PC10", "NK", "Mono", "Gran", "Bcell", "CD8T",
                   "CD4T") # create list of variables to remove NAs


data <- data_unfiltered[complete.cases(data_unfiltered[ , variables_JHS]), ]
dim(data)
# [1] 1227   43

# PLOT RESIDUALS
model_resid <- lm(data$MRS_bonferroni ~ sex + age + currentSmoker + BMI +
                    Income + edu3cat + alc + marital2cat + nbSESanascore +
                    PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + NK + Mono + Gran + 
                    Bcell + CD8T + CD4T, data=data)

data$MRS_resid <- resid(model_resid)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
#1 - 614    2 - 613
v2_bonf_all <- data

cox_model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat, 
                   data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

v2_bonf_cox <- coefficients

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")
write.csv(v2_bonf_cox, file = "v2_cox_results_bonf.csv")

######## MRS_TNF #########
library(dplyr)
library(tidyverse)
library(survival)
library(survminer)

setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/WHI_TNF_sites_JHS")

chd_2 <- read.csv("Input/chd_data_v2.csv", header = TRUE)
mrs_2 <- read.csv("Input/TNF_sites_MRS_JHSv2.csv", header = TRUE)

colnames(mrs_2)[1] <- "SUBJECT_ID"
data_unfiltered <- merge(mrs_2, chd_2, by = "SUBJECT_ID")

data <- data_unfiltered %>% drop_na(CHD) %>% drop_na(days_to_event) %>%
  drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>% drop_na(TNF_MRS)
dim(data)
#[1] 1202   45

res_fit <- lm(TNF_MRS ~ sex + age + 
                currentSmoker + BMI + Income + edu3cat + alc + 
                marital2cat + nbSESanascore + 
                PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                PC8 + PC9 + PC10 + NK + Mono + Gran + 
                Bcell + CD8T + CD4T, data=data)

data$MRS_resid <- resid(res_fit)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat) # 1 - 601,   2 - 601

model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat, data = data)

test_ph <- cox.zph(model)
print(test_ph)

# Summarize the model
coefficients <- summary(model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

write.csv(coefficients, file = "Output/TNF_MRS_CHD_v2.csv")


