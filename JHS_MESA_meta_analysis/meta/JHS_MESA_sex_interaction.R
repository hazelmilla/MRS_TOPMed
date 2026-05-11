# MRS-CHD Sex interaction (Cox Regression)
library(tidyverse)
library(dplyr)
library(survival)

######### JHS v1 ##########

setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/V1/Input_files")
data_unfiltered <- read.csv("chd_data.csv", header = TRUE)
dim(data_unfiltered)
# [1] 1554   40

# Model - sex, age, smoking, bmi, income, education, alc, marital status,
# SES, PC1-10, Houseman cell types
data <- data_unfiltered %>% drop_na(stress2cat) %>% drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>%
  drop_na(PC1) %>% drop_na(PC2) %>% drop_na(PC3) %>% 
  drop_na(PC4) %>% drop_na(PC5) %>% drop_na(PC6) %>% drop_na(PC7) %>% 
  drop_na(PC8) %>% drop_na(PC9) %>% drop_na(PC10) %>% drop_na(NK) %>% 
  drop_na(Mono) %>% drop_na(Gran) %>% drop_na(Bcell) %>% drop_na(CD8T) %>%
  drop_na(CD4T) %>% drop_na(MRS)
dim(data)
# [1] 1312   40

# PLOT RESIDUALS
model_resid <- lm(data$MRS ~ age + currentSmoker + BMI +
                    Income + edu3cat + alc + marital2cat + nbSESanascore +
                    PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + NK + Mono + Gran + 
                    Bcell + CD8T + CD4T, data=data)

data$MRS_resid <- resid(model_resid)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
#1 - 656    2 - 656
JHSv1_FDR_all <- data

cox_model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat*sex, 
                   data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")
write.csv(coefficients, file = "MRS_sex_interaction_CHD_Cox_FDR_v1.csv")


######### JHS v2 ##########

setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/V2")

data_unfiltered <- read.csv("chd_data_v2.csv", header = TRUE)
dim(data_unfiltered)
# [1] 1499   43

# Model - sex, age, smoking, bmi, income, education, alc, marital status,
# SES, PC1-10, Houseman cell types
data <- data_unfiltered %>% drop_na(stress2cat) %>% drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>%
  drop_na(PC1) %>% drop_na(PC2) %>% drop_na(PC3) %>% 
  drop_na(PC4) %>% drop_na(PC5) %>% drop_na(PC6) %>% drop_na(PC7) %>% 
  drop_na(PC8) %>% drop_na(PC9) %>% drop_na(PC10) %>% drop_na(NK) %>% 
  drop_na(Mono) %>% drop_na(Gran) %>% drop_na(Bcell) %>% drop_na(CD8T) %>%
  drop_na(CD4T) %>% drop_na(MRS)
dim(data)
# [1] 1217   43

# PLOT RESIDUALS
model_resid <- lm(data$MRS ~ age + currentSmoker + BMI +
                    Income + edu3cat + alc + marital2cat + nbSESanascore +
                    PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + NK + Mono + Gran + 
                    Bcell + CD8T + CD4T, data=data)

data$MRS_resid <- resid(model_resid)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
#1 - 609    2 - 608
JHSv2_FDR_all <- data

cox_model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat*sex, 
                   data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")
write.csv(coefficients, file = "MRS_sex_interaction_CHD_Cox_FDR_v2.csv")



######### MESA ##########
setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_MESA/Input")
data_unfiltered <- read.csv("prepared_data_MESA.csv", header = TRUE)

var <- c("stress2cat")

cov <- c("AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
         "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
         "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
         "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
         "Bcell", "CD8T", "CD4T")

length(cov)
#28

df <- data_unfiltered[,colnames(data_unfiltered) %in% cov]
dim(df)
#870 28

data <- data_unfiltered[complete.cases(data_unfiltered[ , cov]), ]
data <- data %>% drop_na(stress2cat)
dim(data) #[1] 698  94

resid <- paste("MRS ~", paste(cov, collapse = " + "))
resid_formula <- as.formula(resid)

res_fit <- lm(resid_formula, data=data)

data$MRS_resid <- resid(res_fit)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
# 1 - 349, 2 - 349
mesa_fdr <- data

cox_model <- coxph(Surv(chdatt, chda) ~ MRSresid2cat*gender1, data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

mesa_cox_fdr <- coefficients

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")
write.csv(mesa_cox_fdr, file = "MRS_sex_interaction_CHD_Cox_FDR_MESA.csv")


######### Meta-analysis ##########

library(rmeta)
library(tidyverse)

# JHS data
setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")

v1 <- read.csv("MRS_sex_interaction_CHD_Cox_FDR_v1.csv", header = TRUE)
v2 <- read.csv("MRS_sex_interaction_CHD_Cox_FDR_v2.csv", header = TRUE)
mesa <- read.csv("MRS_sex_interaction_CHD_Cox_FDR_MESA.csv", header = TRUE)

v1$study <- ("EPICv1")
v2$study <- ("EPICv2")
mesa$study <- ("MESA")

meta_cox <- rbind(v1, v2, mesa)
meta_cox <- as.data.frame(meta_cox)

meta_cox$X[9] <- "MRSresid2cat:sexMale"

meta <- meta.summaries(log(HR), SE, method = c("fixed"), logscale = FALSE,
                       names = study, conf.level = 0.95, data = meta_cox,
                       subset = NULL)
run_coef <- c(meta$summary, meta$se.summary, meta$het[3], meta$test, NA)
names(run_coef) <- c("Coefficient", "SE", "p_het", "z", "p_meta", "HR")
run_coef[6] <- exp(run_coef[1])
print(run_coef)

meta_FDR_sex_interaction <- run_coef

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Output")
write.csv(meta_FDR_sex_interaction, "MRS_CHD_sex_interaction_FE_model_meta_JHS_MESA_FDR.csv")




######### JHS v1 ##########

setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/V1/Input_files")
data_unfiltered <- read.csv("chd_data.csv", header = TRUE)
dim(data_unfiltered)
# [1] 1554   40

# Model - sex, age, smoking, bmi, income, education, alc, marital status,
# SES, PC1-10, Houseman cell types
data <- data_unfiltered %>% drop_na(stress2cat) %>% drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>%
  drop_na(PC1) %>% drop_na(PC2) %>% drop_na(PC3) %>% 
  drop_na(PC4) %>% drop_na(PC5) %>% drop_na(PC6) %>% drop_na(PC7) %>% 
  drop_na(PC8) %>% drop_na(PC9) %>% drop_na(PC10) %>% drop_na(NK) %>% 
  drop_na(Mono) %>% drop_na(Gran) %>% drop_na(Bcell) %>% drop_na(CD8T) %>%
  drop_na(CD4T) %>% drop_na(MRS_bonferroni)
dim(data)
# [1] 1312   40

# PLOT RESIDUALS
model_resid <- lm(data$MRS_bonferroni ~ age + currentSmoker + BMI +
                    Income + edu3cat + alc + marital2cat + nbSESanascore +
                    PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + NK + Mono + Gran + 
                    Bcell + CD8T + CD4T, data=data)

data$MRS_resid <- resid(model_resid)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
#1 - 656    2 - 656
JHSv1_FDR_all <- data

cox_model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat*sex, 
                   data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")
write.csv(coefficients, file = "MRS_sex_interaction_CHD_Cox_Bonf_v1.csv")


######### JHS v2 ##########

setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/V2")

data_unfiltered <- read.csv("chd_data_v2.csv", header = TRUE)
dim(data_unfiltered)
# [1] 1499   43

# Model - sex, age, smoking, bmi, income, education, alc, marital status,
# SES, PC1-10, Houseman cell types
data <- data_unfiltered %>% drop_na(stress2cat) %>% drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>%
  drop_na(PC1) %>% drop_na(PC2) %>% drop_na(PC3) %>% 
  drop_na(PC4) %>% drop_na(PC5) %>% drop_na(PC6) %>% drop_na(PC7) %>% 
  drop_na(PC8) %>% drop_na(PC9) %>% drop_na(PC10) %>% drop_na(NK) %>% 
  drop_na(Mono) %>% drop_na(Gran) %>% drop_na(Bcell) %>% drop_na(CD8T) %>%
  drop_na(CD4T) %>% drop_na(MRS_bonferroni)
dim(data)
# [1] 1217   43

# PLOT RESIDUALS
model_resid <- lm(data$MRS_bonferroni ~ age + currentSmoker + BMI +
                    Income + edu3cat + alc + marital2cat + nbSESanascore +
                    PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + NK + Mono + Gran + 
                    Bcell + CD8T + CD4T, data=data)

data$MRS_resid <- resid(model_resid)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
#1 - 609    2 - 608

cox_model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat*sex, 
                   data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")
write.csv(coefficients, file = "MRS_sex_interaction_CHD_Cox_Bonf_v2.csv")



######### MESA ##########
setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_MESA/Input")
data_unfiltered <- read.csv("prepared_data_MESA.csv", header = TRUE)

var <- c("stress2cat")

cov <- c("AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
         "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
         "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
         "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
         "Bcell", "CD8T", "CD4T")

length(cov)
#28

df <- data_unfiltered[,colnames(data_unfiltered) %in% cov]
dim(df)
#870 28

data <- data_unfiltered[complete.cases(data_unfiltered[ , cov]), ]
data <- data %>% drop_na(stress2cat)
dim(data) #[1] 698  94

resid <- paste("MRS_bonf ~", paste(cov, collapse = " + "))
resid_formula <- as.formula(resid)

res_fit <- lm(resid_formula, data=data)

data$MRS_resid <- resid(res_fit)

MRS_median <- summary(data$MRS_resid)[3]
data$MRSresid2cat[data$MRS_resid <= MRS_median] <- 1
data$MRSresid2cat[data$MRS_resid > MRS_median] <- 2
table(data$MRSresid2cat)
# 1 - 349, 2 - 349
mesa_fdr <- data

cox_model <- coxph(Surv(chdatt, chda) ~ MRSresid2cat*gender1, data = data)

test_ph <- cox.zph(cox_model)
print(test_ph)

coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "HR", "SE",
                            "z", "P_val")

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")
write.csv(mesa_cox_fdr, file = "MRS_sex_interaction_CHD_Cox_Bonf_MESA.csv")


######### Meta-analysis ##########

library(rmeta)
library(tidyverse)

# JHS data
setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Input")

v1 <- read.csv("MRS_sex_interaction_CHD_Cox_Bonf_v1.csv", header = TRUE)
v2 <- read.csv("MRS_sex_interaction_CHD_Cox_Bonf_v2.csv", header = TRUE)
mesa <- read.csv("MRS_sex_interaction_CHD_Cox_Bonf_MESA.csv", header = TRUE)

v1$study <- ("EPICv1")
v2$study <- ("EPICv2")
mesa$study <- ("MESA")

meta_cox <- rbind(v1, v2, mesa)
meta_cox <- as.data.frame(meta_cox)

meta_cox$X[9] <- "MRSresid2cat:sexMale"

meta <- meta.summaries(log(HR), SE, method = c("fixed"), logscale = FALSE,
                       names = study, conf.level = 0.95, data = meta_cox,
                       subset = NULL)
run_coef <- c(meta$summary, meta$se.summary, meta$het[3], meta$test, NA)
names(run_coef) <- c("Coefficient", "SE", "p_het", "z", "p_meta", "HR")
run_coef[6] <- exp(run_coef[1])
print(run_coef)

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Output")
write.csv(run_coef, "MRS_CHD_sex_interaction_FE_model_meta_JHS_MESA_Bonferroni.csv")
