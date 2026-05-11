############## MRS permutations ################

library(dplyr)
library(tidyverse)

###### JHS V1 ######

# start with FDR CpGs
setwd("/work/users/h/e/hemilla/JHS_stress_CpG/Input_files")
meta_coef <- read.csv("meta_coef_sig_annotated_SLEs.csv", header = TRUE)
head(meta_coef)

load("beta.noXY.RData")
load("beta.v2.subjid.RData")

setwd("/work/users/h/e/hemilla/WHI_MESA_MRS/Input")
load("MESA_methyl_visit_1_cleaned.RData")

# get only CpG beta values
beta.mesa <- methylation.vis1[, grep("^cg", colnames(methylation.vis1))]
beta.mesa <- t(beta.mesa) # transpose beta so rows are CpGs and cols are Sample ID

meta_coef_bonf <- meta_coef[meta_coef$bonferroni < 0.05,] 

MRS_weights <- meta_coef$beta # take MRS weights for FDR sites
MRS_weights_bonf <- meta_coef_bonf$beta # take MRS weights for Bonferroni sites

names(MRS_weights) <- meta_coef$IlmnID
names(MRS_weights_bonf) <- meta_coef_bonf$IlmnID

# This will be all sites for cohorts, not just WHI FDR significant sites
common_probes <- Reduce(intersect, list(rownames(beta.noXY), rownames(beta), rownames(beta.mesa)))
length(common_probes) #[1] 659563

set.seed(123)
all_cpgs <- common_probes

n_permutations <- 10000
random_cpg_sets <- replicate(n_permutations, sample(all_cpgs, 841), simplify=FALSE)
random_cpg_sets_bonf <- replicate(n_permutations, sample(all_cpgs, 13), simplify=FALSE)

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
saveRDS(random_cpg_sets, "random_cpg_sets.rds")
saveRDS(random_cpg_sets_bonf, "random_cpg_sets_bonf.rds")

######## Permutation calculation ##########

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
random_cpg_sets <- readRDS("random_cpg_sets.rds")

# JHS V1
beta.v1 <- beta.noXY[rownames(beta.noXY) %in% common_probes,]
beta.v1 <- t(beta.v1)

null_scores_matrix <- matrix(NA, nrow=nrow(beta.v1), ncol = n_permutations)
rownames(null_scores_matrix) <- rownames(beta.v1)

for (i in 1:n_permutations){
  random_cpgs <- random_cpg_sets[[i]]
  random_betas <- beta.v1[,random_cpgs]
  
  null_score <- as.vector(as.matrix(random_betas) %*% MRS_weights)
  null_scores_matrix[,i] <- null_score
}

saveRDS(null_scores_matrix, file = "JHSv1_MRS_permutations_FDR.rds")

# JHS V2
beta.v2 <- beta[rownames(beta) %in% common_probes,]
beta.v2 <- t(beta.v2)

null_scores_matrix <- matrix(NA, nrow=nrow(beta.v2), ncol = n_permutations)
rownames(null_scores_matrix) <- rownames(beta.v2)

for (i in 1:n_permutations){
  random_cpgs <- random_cpg_sets[[i]]
  random_betas <- beta.v2[,random_cpgs]
  
  null_score <- as.vector(as.matrix(random_betas) %*% MRS_weights)
  null_scores_matrix[,i] <- null_score
}

saveRDS(null_scores_matrix, file = "JHSv2_MRS_permutations_FDR.rds")

# MESA
beta.mesa1 <- beta.mesa[rownames(beta.mesa) %in% common_probes,]
beta.mesa1 <- t(beta.mesa1)

null_scores_matrix <- matrix(NA, nrow=nrow(beta.mesa1), ncol = n_permutations)
rownames(null_scores_matrix) <- rownames(beta.mesa1)

for (i in 1:n_permutations){
  random_cpgs <- random_cpg_sets[[i]]
  random_betas <- beta.mesa1[,random_cpgs]
  
  null_score <- as.vector(as.matrix(random_betas) %*% MRS_weights)
  null_scores_matrix[,i] <- null_score
}

saveRDS(null_scores_matrix, file = "MESA_MRS_permutations_FDR.rds")

###### Residualize MRS #########

library(dplyr)
library(tidyverse)

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
null_scores_matrix <- readRDS("JHSv1_MRS_permutations_FDR.rds")

setwd("/work/users/h/e/hemilla/JHS_stress_CpG/Input_files")
covariates <- read.csv("prepared_data.csv", header = TRUE)

covariates <- covariates %>% drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>%
  drop_na(PC1) %>% drop_na(PC2) %>% drop_na(PC3) %>% 
  drop_na(PC4) %>% drop_na(PC5) %>% drop_na(PC6) %>% drop_na(PC7) %>% 
  drop_na(PC8) %>% drop_na(PC9) %>% drop_na(PC10) %>% drop_na(NK) %>% 
  drop_na(Mono) %>% drop_na(Gran) %>% drop_na(Bcell) %>% drop_na(CD8T) %>%
  drop_na(CD4T)
dim(covariates) #1427

rownames(covariates) <- covariates$Sample_Name
common_ids <- intersect(rownames(null_scores_matrix), rownames(covariates))
null_scores_matrix <- null_scores_matrix[common_ids, , drop = FALSE]
null_scores_matrix <- null_scores_matrix[rownames(null_scores_matrix) %in% covariates$Sample_Name, ,drop = FALSE]

null_scores_matrix <- null_scores_matrix[order(match(rownames(null_scores_matrix), covariates$Sample_Name)),]
identical((rownames(null_scores_matrix) %>% as.integer()), covariates$Sample_Name)
dim(null_scores_matrix) #[1]    1427 10000

resid_matrix <- matrix(NA, nrow = nrow(null_scores_matrix), ncol = ncol(null_scores_matrix))

for (i in 1:ncol(null_scores_matrix)){
  # Extract permuted MRS
  perm_score <- null_scores_matrix[,i]
  covariates$perm_score <- perm_score
  
  # Fit linear model to get residuals
  fit <- lm(perm_score ~ sex + age + 
                currentSmoker + BMI + Income + edu3cat + alc + 
                marital2cat + nbSESanascore + 
                PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                PC8 + PC9 + PC10 + NK + Mono + Gran + 
                Bcell + CD8T + CD4T, data = covariates)
  
  # Store residuals 
  resid_matrix[,i] <- resid(fit)
}

rownames(resid_matrix) <- covariates$Sample_Name

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
saveRDS(resid_matrix, file = "JHSv1_MRS_permutations_residualized_FDR.rds")

# JHS V2
setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
null_scores_matrix <- readRDS("JHSv2_MRS_permutations_FDR.rds")

setwd("/work/users/h/e/hemilla/JHS_stress_CpG/Input_files")
covariates <- read.csv("prepared_data_v2.csv", header = TRUE)

covariates <- covariates %>% drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>%
  drop_na(PC1) %>% drop_na(PC2) %>% drop_na(PC3) %>% 
  drop_na(PC4) %>% drop_na(PC5) %>% drop_na(PC6) %>% drop_na(PC7) %>% 
  drop_na(PC8) %>% drop_na(PC9) %>% drop_na(PC10) %>% drop_na(NK) %>% 
  drop_na(Mono) %>% drop_na(Gran) %>% drop_na(Bcell) %>% drop_na(CD8T) %>%
  drop_na(CD4T)
dim(covariates) #[1] 1323   36

rownames(covariates) <- covariates$SUBJECT_ID
common_ids <- intersect(rownames(null_scores_matrix), rownames(covariates))
null_scores_matrix <- null_scores_matrix[common_ids, , drop = FALSE]
null_scores_matrix <- null_scores_matrix[rownames(null_scores_matrix) %in% covariates$SUBJECT_ID, ,drop = FALSE]

null_scores_matrix <- null_scores_matrix[order(match(rownames(null_scores_matrix), covariates$SUBJECT_ID)),]
identical((rownames(null_scores_matrix) %>% as.integer()), covariates$SUBJECT_ID)
dim(null_scores_matrix) #[1]    1323 10000

resid_matrix <- matrix(NA, nrow = nrow(null_scores_matrix), ncol = ncol(null_scores_matrix))

for (i in 1:ncol(null_scores_matrix)){
  # Extract permuted MRS
  perm_score <- null_scores_matrix[,i]
  covariates$perm_score <- perm_score
  
  # Fit linear model to get residuals
  fit <- lm(perm_score ~ sex + age + 
              currentSmoker + BMI + Income + edu3cat + alc + 
              marital2cat + nbSESanascore + 
              PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + NK + Mono + Gran + 
              Bcell + CD8T + CD4T, data = covariates)
  
  # Store residuals 
  resid_matrix[,i] <- resid(fit)
}

rownames(resid_matrix) <- covariates$SUBJECT_ID

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
saveRDS(resid_matrix, file = "JHSv2_MRS_permutations_residualized_FDR.rds")

# MESA
setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
null_scores_matrix <- readRDS("MESA_MRS_permutations_FDR.rds")

setwd("/work/users/h/e/hemilla/WHI_MESA_MRS/Input")
covariates <- read.csv("prepared_data_MESA.csv", header = TRUE)

cov <- c("gender1", "AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
         "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
         "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
         "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
         "Bcell", "CD8T", "CD4T")
covariates <- covariates[complete.cases(covariates[ , cov]), ]
dim(covariates) #[1] 610  93

rownames(covariates) <- covariates$TOEID
common_ids <- intersect(rownames(null_scores_matrix), rownames(covariates))
null_scores_matrix <- null_scores_matrix[common_ids, , drop = FALSE]
null_scores_matrix <- null_scores_matrix[rownames(null_scores_matrix) %in% covariates$TOEID, ,drop = FALSE]

null_scores_matrix <- null_scores_matrix[order(match(rownames(null_scores_matrix), covariates$TOEID)),]
identical(rownames(null_scores_matrix), covariates$TOEID)
dim(null_scores_matrix) #[1]    610 10000

resid_matrix <- matrix(NA, nrow = nrow(null_scores_matrix), ncol = ncol(null_scores_matrix))

for (i in 1:ncol(null_scores_matrix)){
  # Extract permuted MRS
  perm_score <- null_scores_matrix[,i]
  covariates$perm_score <- perm_score
  
  # Calculate residuals
  resid_string <- paste("perm_score ~", paste(cov, collapse = " + "))
  resid_formula <- as.formula(resid_string)
  
  fit <- lm(resid_formula, data=covariates)
  
  # Store residuals 
  resid_matrix[,i] <- resid(fit)
}

rownames(resid_matrix) <- covariates$TOEID

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
saveRDS(resid_matrix, file = "MESA_MRS_permutations_residualized_FDR.rds")

#### Determine association with CHD ####
library(survival)

# JHS V1
setwd("/work/users/h/e/hemilla/JHS_stress_CpG/Input_files")
covariates <- read.csv("chd_data.csv", header = TRUE)

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
perm_matrix <- readRDS("JHSv1_MRS_permutations_residualized_FDR.rds")

covariates <- covariates %>% drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>%
  drop_na(PC1) %>% drop_na(PC2) %>% drop_na(PC3) %>% 
  drop_na(PC4) %>% drop_na(PC5) %>% drop_na(PC6) %>% drop_na(PC7) %>% 
  drop_na(PC8) %>% drop_na(PC9) %>% drop_na(PC10) %>% drop_na(NK) %>% 
  drop_na(Mono) %>% drop_na(Gran) %>% drop_na(Bcell) %>% drop_na(CD8T) %>%
  drop_na(CD4T)
dim(covariates) #1321

rownames(perm_matrix)
rownames(covariates) <- covariates$Sample_Name

common_ids <- intersect(rownames(perm_matrix), rownames(covariates))
perm_matrix <- perm_matrix[common_ids, , drop = FALSE]
perm_matrix <- perm_matrix[rownames(perm_matrix) %in% covariates$Sample_Name, ,drop = FALSE]

perm_matrix <- perm_matrix[order(match(rownames(perm_matrix), covariates$Sample_Name)),]
identical((rownames(perm_matrix) %>% as.integer()), covariates$Sample_Name)
dim(perm_matrix) #[1]    1321 10000

n_perm <- ncol(perm_matrix)
perm_results <- matrix(NA, nrow = n_perm, ncol = 5)
colnames(perm_results) <- c("Iteration", "HR", "SE", "Z", "P_value")
ph_tests <- list()

for (i in 1:n_perm) {
  covariates$MRS_perm <- perm_matrix[,i]
  
  MRS_median <- summary(covariates$MRS_perm)[3]
  covariates$MRSresid2cat[covariates$MRS_perm <= MRS_median] <- 1
  covariates$MRSresid2cat[covariates$MRS_perm > MRS_median] <- 2
  table(covariates$MRSresid2cat)
  
  cox_model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat, data = covariates)
  
  coef_summary <- summary(cox_model)$coefficients
  
  perm_results[i,] <- c(i,
                        coef_summary["MRSresid2cat", "exp(coef)"],
                        coef_summary["MRSresid2cat", "se(coef)"],
                        coef_summary["MRSresid2cat", "z"],
                        coef_summary["MRSresid2cat", "Pr(>|z|)"])
  
  ph_tests[[i]] <- cox.zph(cox_model)
  
}  

perm_results_df <- perm_results %>% as.data.frame()

#check HR value ranges
range(perm_results_df$HR) #[1] 0.6208869 1.5278398
sum(perm_results_df$P_value < 0.05) #13

names(ph_tests) <- 1:1000
violating_vars <- names(ph_tests)[sapply(ph_tests, function(x) x$table[1, 3] < 0.05)]
print(violating_vars)

write.csv(perm_results_df, file = "JHSv1_permutation_results.csv", row.names = FALSE)

# JHS V2
setwd("/work/users/h/e/hemilla/JHS_stress_CpG/Input_files")
covariates <- read.csv("chd_data_v2.csv", header = TRUE)

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
perm_matrix <- readRDS("JHSv2_MRS_permutations_residualized_FDR.rds")

covariates <- covariates %>% drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>%
  drop_na(PC1) %>% drop_na(PC2) %>% drop_na(PC3) %>% 
  drop_na(PC4) %>% drop_na(PC5) %>% drop_na(PC6) %>% drop_na(PC7) %>% 
  drop_na(PC8) %>% drop_na(PC9) %>% drop_na(PC10) %>% drop_na(NK) %>% 
  drop_na(Mono) %>% drop_na(Gran) %>% drop_na(Bcell) %>% drop_na(CD8T) %>%
  drop_na(CD4T)
dim(covariates) #1227

rownames(perm_matrix)
rownames(covariates) <- covariates$SUBJECT_ID

common_ids <- intersect(rownames(perm_matrix), rownames(covariates))
perm_matrix <- perm_matrix[common_ids, , drop = FALSE]
perm_matrix <- perm_matrix[rownames(perm_matrix) %in% covariates$SUBJECT_ID, ,drop = FALSE]

perm_matrix <- perm_matrix[order(match(rownames(perm_matrix), covariates$SUBJECT_ID)),]
identical((rownames(perm_matrix) %>% as.integer()), covariates$SUBJECT_ID)
dim(perm_matrix) #[1] 1227 10000

n_perm <- ncol(perm_matrix)
perm_results <- matrix(NA, nrow = n_perm, ncol = 5)
colnames(perm_results) <- c("Iteration", "HR", "SE", "Z", "P_value")
ph_tests <- list()

for (i in 1:n_perm) {
  covariates$MRS_perm <- perm_matrix[,i]
  
  MRS_median <- summary(covariates$MRS_perm)[3]
  covariates$MRSresid2cat[covariates$MRS_perm <= MRS_median] <- 1
  covariates$MRSresid2cat[covariates$MRS_perm > MRS_median] <- 2
  table(covariates$MRSresid2cat)
  
  cox_model <- coxph(Surv(days_to_event, CHD) ~ MRSresid2cat, data = covariates)
  
  coef_summary <- summary(cox_model)$coefficients
  
  perm_results[i,] <- c(i,
                        coef_summary["MRSresid2cat", "exp(coef)"],
                        coef_summary["MRSresid2cat", "se(coef)"],
                        coef_summary["MRSresid2cat", "z"],
                        coef_summary["MRSresid2cat", "Pr(>|z|)"])
  
  ph_tests[[i]] <- cox.zph(cox_model)
  
}  

perm_results_df <- perm_results %>% as.data.frame()

#check HR value ranges
range(perm_results_df$HR)
sum(perm_results_df$P_value < 0.05) #20

names(ph_tests) <- 1:1000
violating_vars <- names(ph_tests)[sapply(ph_tests, function(x) x$table[1, 3] < 0.05)]
print(violating_vars)

write.csv(perm_results_df, file = "JHSv2_permutation_results.csv", row.names = FALSE)

# MESA
setwd("/work/users/h/e/hemilla/WHI_MESA_MRS/Input")
covariates <- read.csv("prepared_data_MESA.csv", header = TRUE)

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
perm_matrix <- readRDS("MESA_MRS_permutations_residualized_FDR.rds")

cov <- c("gender1", "AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
         "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
         "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
         "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
         "Bcell", "CD8T", "CD4T", "chda")
covariates <- covariates[complete.cases(covariates[ , cov]), ]
dim(covariates) #[1] 609  93

rownames(perm_matrix)
rownames(covariates) <- covariates$TOEID

common_ids <- intersect(rownames(perm_matrix), rownames(covariates))
perm_matrix <- perm_matrix[common_ids, , drop = FALSE]
perm_matrix <- perm_matrix[rownames(perm_matrix) %in% covariates$TOEID, ,drop = FALSE]

perm_matrix <- perm_matrix[order(match(rownames(perm_matrix), covariates$TOEID)),]
identical(rownames(perm_matrix), covariates$TOEID)
dim(perm_matrix) #[1] 707

n_perm <- ncol(perm_matrix)
perm_results <- matrix(NA, nrow = n_perm, ncol = 5)
colnames(perm_results) <- c("Iteration", "HR", "SE", "Z", "P_value")
ph_tests <- list()

for (i in 1:n_perm) {
  covariates$MRS_perm <- perm_matrix[,i]
  
  MRS_median <- summary(covariates$MRS_perm)[3]
  covariates$MRSresid2cat[covariates$MRS_perm <= MRS_median] <- 1
  covariates$MRSresid2cat[covariates$MRS_perm > MRS_median] <- 2
  table(covariates$MRSresid2cat)
  
  cox_model <- coxph(Surv(chdatt, chda) ~ MRSresid2cat, data = covariates)
  
  coef_summary <- summary(cox_model)$coefficients
  
  perm_results[i,] <- c(i,
                        coef_summary["MRSresid2cat", "exp(coef)"],
                        coef_summary["MRSresid2cat", "se(coef)"],
                        coef_summary["MRSresid2cat", "z"],
                        coef_summary["MRSresid2cat", "Pr(>|z|)"])
  
  ph_tests[[i]] <- cox.zph(cox_model)
  
}  

perm_results_df <- perm_results %>% as.data.frame()

#check HR value ranges
range(perm_results_df$HR) #0.6032669 2.6788420
sum(perm_results_df$P_value < 0.05) #314

names(ph_tests) <- 1:1000
violating_vars <- names(ph_tests)[sapply(ph_tests, function(x) x$table[1, 3] < 0.05)]
print(violating_vars)

write.csv(perm_results_df, file = "MESA_permutation_results.csv", row.names = FALSE)

####### Permutation meta-analysis ########

library(rmeta)
library(tidyverse)
library(dplyr)

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")

v1 <- read.csv("JHSv1_permutation_results.csv", header = TRUE)
v2 <- read.csv("JHSv2_permutation_results.csv", header = TRUE)
mesa <- read.csv("MESA_permutation_results.csv", header = TRUE)

v1$study <- ("EPICv1")
v2$study <- ("EPICv2")
mesa$study <- ("MESA")

sitebind <- NULL
run_meta <- NULL
run_coef <- NULL
meta_coef <- NULL
i <- 1

for (i in (1:length(v1$Iteration))){
  meta_cox <- rbind(v1[i,], v2[i,], mesa[i,])
  meta <- meta.summaries(log(HR), SE, method = c("fixed"), logscale = FALSE,
                         names = study, conf.level = 0.95, data = meta_cox,
                         subset = NULL)
  run_coef <- c(meta$summary, meta$se.summary, meta$het[3], meta$test, NA)
  meta_coef <- rbind(meta_coef, run_coef)
  print(i)
  flush.console()
}

dim(meta_coef)
class(meta_coef)
meta_coef_df <- as.data.frame(meta_coef)
colnames(meta_coef_df) <- c("Coefficient", "SE", "p_het", "z", "p_meta", "HR")
rownames(meta_coef_df) <- v1$CpG

meta_coef_df$HR <- exp(meta_coef_df$Coefficient)

meta_coef_df$q_meta <- p.adjust(meta_coef_df$p_meta, method = "fdr")
meta_coef_df$bonferroni <- p.adjust(meta_coef_df$p_meta, method = "bonferroni")

sum(meta_coef_df$p_het < 0.05) # 98
sum(meta_coef_df$p_meta < 0.05) # 9
sum(meta_coef_df$q_meta < 0.05) # 0
sum(meta_coef_df$bonferroni < 0.05) # 0

write.csv(meta_coef_df, file = "MRS_permutations_CHD_MESA_JHS_meta.csv")

meta_coef_df_non_hg <- meta_coef_df %>% filter(p_het > 0.05)
sum(meta_coef_df_non_hg$p_meta < 0.05) # 1

### Heterogeneous permutations ###
sig_hg <- meta_coef_df %>% filter(meta_coef_df$p_het < 0.05)
sig_hg <- rownames(sig_hg)

mesa_hg <- mesa[mesa$Iteration %in% sig_hg,]
v1_hg <- v1[v1$Iteration %in% sig_hg,]
v2_hg <- v2[v2$Iteration %in% sig_hg,]

sitebind <- NULL
run_meta <- NULL
run_coef <- NULL
meta_coef <- NULL
i <- 1

for (i in (1:length(v1_hg$Iteration))){
  meta_cox <- rbind(v1_hg[i,], v2_hg[i,], mesa_hg[i,])
  meta <- meta.summaries(log(HR), SE, method = c("random"), logscale = FALSE,
                         names = study, conf.level = 0.95, data = meta_cox,
                         subset = NULL)
  run_coef <- c(meta$summary, meta$se.summary, meta$het[3], meta$test, NA)
  meta_coef <- rbind(meta_coef, run_coef)
  print(i)
  flush.console()
}

dim(meta_coef)
class(meta_coef)
meta_coef_df <- as.data.frame(meta_coef)
colnames(meta_coef_df) <- c("Coefficient", "SE", "p_het", "z", "p_meta", "HR")
rownames(meta_coef_df) <- v1_hg$Iteration

meta_coef_df$HR <- exp(meta_coef_df$Coefficient)

meta_coef_df$q_meta <- p.adjust(meta_coef_df$p_meta, method = "fdr")
meta_coef_df$bonferroni <- p.adjust(meta_coef_df$p_meta, method = "bonferroni")

sum(meta_coef_df$p_het < 0.05) # 1
sum(meta_coef_df$p_meta < 0.05) # 0
sum(meta_coef_df$q_meta < 0.05) # 0
sum(meta_coef_df$bonferroni < 0.05) # 0

write.csv(meta_coef_df, file = "MRS_permutations_CHD_MESA_JHS_heterogeneous_iterations_meta.csv")


# Calculate p-value
setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
meta_coef_df <- read.csv("MRS_permutations_CHD_MESA_JHS_meta.csv")
meta_coef_df_non_hg <- read.csv("MRS_permutations_CHD_MESA_JHS_heterogeneous_iterations_meta.csv")

z_obs <- 1.1225
z_perm <- c(meta_coef_df$z, meta_coef_df_non_hg$z)
p_perm <- (sum(z_perm >= z_obs) + 1) / (length(z_perm) + 1) #[1] 0.01703139

library(ggplot2)
df <- data.frame(z = z_perm)

z_plot <- ggplot(df, aes(x = z)) + 
  geom_density(fill = "lightgrey", alpha = 0.6, color = "darkgrey", linewidth = 1.2) +
  geom_vline(xintercept = z_obs, color = "#d73027", linetype = "dashed", linewidth = 1) +
  annotate("text", 
           x = z_obs,
           y = max(density(z_perm)$y),
           label = paste("Observed z =", round(z_obs, 2), "\n", "p-value =", round(p_perm, 2)),
           hjust = -0.1,
           vjust = 0.5,
           color = "#d73027", 
           size = 4.2) +
  labs(
    title = "Permutation Z-score Distribution",
    x = "z-score",
    y = "Density"
) + 
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

z_plot

ggsave("JHS_MESA_MRS_CHD_permutations_FDR.png", z_plot, width = 10, 
       height = 7, units = "in", dpi = 300)
