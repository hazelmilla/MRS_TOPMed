############## MRS permutation resampling analysis ################

library(dplyr)
library(tidyverse)

###### JHS V1 ######

# start with FDR CpGs
meta_coef <- read.csv("/proj/azannas/projects/Hazel/JHS/Input/meta_coef_sig_annotated_SLEs.csv", header = TRUE)
head(meta_coef)
load("/proj/azannas/projects/Hazel/JHS/Input/beta.v2.subjid.RData")
load("/proj/azannas/projects/DNAm_Data_JHSCC/beta.noXY.RData")
load("/work/users/h/e/hemilla/MESA/Input/MESA_methyl_visit_1_cleaned.RData")

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

setwd("/work/users/h/e/hemilla/WHI_validation/MRS_permutations/")
saveRDS(random_cpg_sets, "random_cpg_sets.rds")
saveRDS(random_cpg_sets_bonf, "random_cpg_sets_bonf.rds")

######## Permutation calculation ##########

setwd("/work/users/h/e/hemilla/WHI_validation/MRS_permutations/")
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

######## MRS Bonferroni Permutations ##########
library(dplyr)
library(tidyverse)

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

n_permutations <- 10000

# Run permutations

setwd("/work/users/h/e/hemilla/MRS_permutations/Output")
random_cpg_sets <- readRDS("random_cpg_sets_bonf.rds")

# JHS V1
beta.v1 <- beta.noXY[rownames(beta.noXY) %in% common_probes,]
beta.v1 <- t(beta.v1)

null_scores_matrix <- matrix(NA, nrow=nrow(beta.v1), ncol = n_permutations)
rownames(null_scores_matrix) <- rownames(beta.v1)

for (i in 1:n_permutations){
  random_cpgs <- random_cpg_sets[[i]]
  random_betas <- beta.v1[,random_cpgs]
  
  null_score <- as.vector(as.matrix(random_betas) %*% MRS_weights_bonf)
  null_scores_matrix[,i] <- null_score
}

saveRDS(null_scores_matrix, file = "JHSv1_MRS_permutations_Bonferroni.rds")

# JHS V2
beta.v2 <- beta[rownames(beta) %in% common_probes,]
beta.v2 <- t(beta.v2)

null_scores_matrix <- matrix(NA, nrow=nrow(beta.v2), ncol = n_permutations)
rownames(null_scores_matrix) <- rownames(beta.v2)

for (i in 1:n_permutations){
  random_cpgs <- random_cpg_sets[[i]]
  random_betas <- beta.v2[,random_cpgs]
  
  null_score <- as.vector(as.matrix(random_betas) %*% MRS_weights_bonf)
  null_scores_matrix[,i] <- null_score
}

saveRDS(null_scores_matrix, file = "JHSv2_MRS_permutations_Bonferroni.rds")

# MESA
beta.mesa1 <- beta.mesa[rownames(beta.mesa) %in% common_probes,]
beta.mesa1 <- t(beta.mesa1)

null_scores_matrix <- matrix(NA, nrow=nrow(beta.mesa1), ncol = n_permutations)
rownames(null_scores_matrix) <- rownames(beta.mesa1)

for (i in 1:n_permutations){
  random_cpgs <- random_cpg_sets[[i]]
  random_betas <- beta.mesa1[,random_cpgs]
  
  null_score <- as.vector(as.matrix(random_betas) %*% MRS_weights_bonf)
  null_scores_matrix[,i] <- null_score
}

saveRDS(null_scores_matrix, file = "MESA_MRS_permutations_Bonferroni.rds")

#######################################################
########### Test null score CHD prediction ############
#######################################################

###### Residualize MRS, run Cox for permuted MRSs, plot z-scores #########

library(dplyr)
library(tidyverse)
library(survival)
library(ggplot2)
library(ggtext)
library(patchwork)
library(rmeta)

setwd("/work/users/h/e/hemilla/WHI_validation/MRS_permutations/")

########################################
############## MRS 841 #################
########################################

####### Prep data #######
null_scores_matrix1 <- readRDS("JHSv1_MRS_permutations_FDR.rds")
null_scores_matrix2 <- readRDS("JHSv2_MRS_permutations_FDR.rds")

jhs_v1 <- read.csv("/proj/azannas/projects/Hazel/TOPMED_MRS/JHS_MRS_data/prepared_data.csv", header = TRUE)
jhs_v2 <- read.csv("/proj/azannas/projects/Hazel/TOPMED_MRS/JHS_MRS_data/prepared_data_v2.csv", header = TRUE)

jhs_cov <- c("sex", "age", "currentSmoker", "BMI", "Income",
             "edu3cat", "alc", "marital2cat", "nbSESanascore",
             "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",
             "PC9", "PC10", "NK", "Mono", "Gran", "Bcell", "CD8T",
             "CD4T")

null_scores_matrix_mesa <- readRDS("MESA_MRS_permutations_FDR.rds")
mesa <- read.csv("/work/users/h/e/hemilla/MESA/Input/prepared_data_MESA.csv", header = TRUE)
mesa_cov <- c("gender1", "AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
              "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
              "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
              "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
              "Bcell", "CD8T", "CD4T")

mesa <- mesa[complete.cases(mesa[ , mesa_cov]), ]
nrow(mesa) # n = 708

jhs_v1 <- jhs_v1[complete.cases(jhs_v1[ , jhs_cov]), ]
nrow(jhs_v1) # n = 1427

jhs_v2 <- jhs_v2[complete.cases(jhs_v2[ , jhs_cov]), ]
nrow(jhs_v2) # n = 1323

######## Residualize MRS #########

get_null_matrix <- function(null_scores_matrix, sample_name, df){
  data <- df
  sample_names <- data[[sample_name]]
  
  print(dim(null_scores_matrix))
  
  # Make sure both are the same type of object
  rownames(null_scores_matrix) <- rownames(null_scores_matrix) %>% as.character()
  sample_names <- as.character(sample_names)
  
  null_scores_matrix <- null_scores_matrix[rownames(null_scores_matrix) %in% sample_names, , drop = FALSE]
  null_scores_matrix <- null_scores_matrix[order(match(rownames(null_scores_matrix), sample_names)),]
  
  print(identical((rownames(null_scores_matrix) %>% as.character()), sample_names))
  print(dim(null_scores_matrix))
  
  return(null_scores_matrix)
}

dataset_names <- c("jhs_v1", "jhs_v2", "mesa")
datasets <- as.list(mget(dataset_names))
names(datasets) <- dataset_names

matrix_names <- c("null_scores_matrix1", "null_scores_matrix2",
                  "null_scores_matrix_mesa")
matrices <- as.list(mget(matrix_names))
names(matrices) <- dataset_names

subjid_names <- list("Sample_Name", "SUBJECT_ID", "TOEID")
names(subjid_names) <- dataset_names

null_matrices <- lapply(dataset_names, function(dataset){
  get_null_matrix(
    null_scores_matrix = matrices[[dataset]],
    sample_name = subjid_names[[dataset]],
    df = datasets[[dataset]]
  )
})

#[1] V1 dim before filtering: 1709 10000
#[1] TRUE
#[1] V1 dim after filtering: 1427 10000
#[1] V2 dim before filtering: 1687 10000
#[1] TRUE
#[1] V2 dim after filtering: 1323 10000
#[1] MESA dim before filtering: 870 10000
#[1] TRUE
#[1] MESA dim after filtering: 708 10000

names(null_matrices) <- matrix_names
list2env(null_matrices, envir = .GlobalEnv)

matrix_names <- c("null_scores_matrix1", "null_scores_matrix2",
                  "null_scores_matrix_mesa")
matrices <- as.list(mget(matrix_names))
names(matrices) <- dataset_names

get_resid <- function(null_scores_matrix, covariates, df){
  # Initialize resid matrix
  resid_matrix <- matrix(NA, nrow = nrow(null_scores_matrix), ncol = ncol(null_scores_matrix))
  rownames(resid_matrix) <- rownames(null_scores_matrix)
  colnames(resid_matrix) <- colnames(null_scores_matrix)
  
  for (i in seq_len(ncol(null_scores_matrix))){
    
    # Extract permuted MRS
    perm_score <- null_scores_matrix[,i]
    data <- df
    data$perm_score <- perm_score
    
    # Fit linear model to get residuals
    formula_text <- paste0("perm_score ~ ", paste(covariates, collapse = " + "))
    formula <- as.formula(formula_text)
    fit <- lm(formula_text, data = data)
    
    # Store residuals 
    resid_matrix[,i] <- resid(fit)
    }
  
  return(resid_matrix)
}

cov_names <- c("jhs_cov", "jhs_cov",
               "mesa_cov")
cov <- mget(cov_names)
names(cov) <- dataset_names

residuals <- lapply(dataset_names, function(dataset){
  get_resid(
  null_scores_matrix = matrices[[dataset]],
  covariates = cov[[dataset]],
  df = datasets[[dataset]]
  )
})

names(residuals) <- paste0(dataset_names, "_resid")
list2env(residuals, envir = .GlobalEnv)

# Save the list
saveRDS(residuals, file = "perm_resid_data_FDR.rds")

#### Determine association with CHD ####

# Load resid data
loaded_list <- readRDS("perm_resid_data_FDR.rds")
dataset_names <- c("jhs_v1", "jhs_v2", "mesa")
names(loaded_list) <- paste0(dataset_names, "_resid")
list2env(loaded_list, envir = .GlobalEnv)

# JHS V1
jhs_v1_chd <- read.csv("/proj/azannas/projects/Hazel/JHS/Input/chd_data.csv", header = TRUE)
jhs_v2_chd <- read.csv("/proj/azannas/projects/Hazel/JHS/Input/chd_data_v2.csv", header = TRUE)
all_jhs_chd_dat <- read.csv("/proj/azannas/projects/JHS/incevtchd_2020_update_dbgap_id.csv")

jhs_v1_chd <- jhs_v1_chd %>%
  left_join(all_jhs_chd_dat, by = "SUBJECT_ID", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y"))

jhs_v2_chd <- jhs_v2_chd %>%
  left_join(all_jhs_chd_dat, by = "SUBJECT_ID", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y"))

dim(jhs_v1_chd) #n = 1554
dim(jhs_v2_chd) #n = 1499

mesa_chd <- filter(mesa, exall != 1)
mesa_chd <- mesa_chd[complete.cases(mesa_chd[ , mesa_cov]), ]
nrow(mesa_chd) # n = 707

jhs_v1_chd <- jhs_v1_chd[complete.cases(jhs_v1_chd[ , jhs_cov]), ]
nrow(jhs_v1_chd) # n = 1321

jhs_v2_chd <- jhs_v2_chd[complete.cases(jhs_v2_chd[ , jhs_cov]), ]
nrow(jhs_v2_chd) # n = 1227

get_perm_resid_dat <- function(resid_matrix, sample_name, df){
  data <- df
  sample_names <- as.character(data[[sample_name]])
  
  common_ids <- intersect(rownames(resid_matrix), sample_names)
  rownames(resid_matrix) <- as.character(rownames(resid_matrix))
  resid_matrix <- resid_matrix[as.character(common_ids), , drop = FALSE]

  resid_matrix <- resid_matrix[order(match(rownames(resid_matrix), sample_names)),]
  print(identical(rownames(resid_matrix), sample_names))
  
  print(dim(resid_matrix))
  
  return(resid_matrix)
}

dataset_names <- c("jhs_v1_chd", "jhs_v2_chd", "mesa_chd")
datasets <- as.list(mget(dataset_names))
names(datasets) <- dataset_names

matrix_names <- names(loaded_list)
matrices <- as.list(mget(matrix_names))
names(matrices) <- dataset_names

subjid_names <- list("Sample_Name", "SUBJECT_ID", "TOEID")
names(subjid_names) <- dataset_names

perm_resid_dat <- lapply(dataset_names, function(dataset){
  get_perm_resid_dat(
    resid_matrix = matrices[[dataset]],
    sample_name = subjid_names[[dataset]],
    df = datasets[[dataset]]
  )
})

#[1] TRUE
#[1] V1 dim after filtering: 1321 10000
#[1] TRUE
#[1] V2 dim after filtering: 1227 10000
#[1] TRUE
#[1] MESA dim after filtering: 707 10000

names(perm_resid_dat) <- matrix_names
list2env(perm_resid_dat, envir = .GlobalEnv)

cox_fun <- function(resid_matrix, df, x, y, time_to){
  # initialize results matrix
  n_perm <- ncol(resid_matrix)
  results <- matrix(NA, nrow = n_perm, ncol = 5)
  colnames(results) <- c("Iteration", "HR", "SE", "Z", "P_value")
  ph_tests <- list()
  
  data <- df
  
  for (i in 1:n_perm) {
    data$MRS_perm <- resid_matrix[,i]
    
    MRS_median <- summary(data$MRS_perm)[3]
    data$MRSresid2cat[data$MRS_perm <= MRS_median] <- 1
    data$MRSresid2cat[data$MRS_perm > MRS_median] <- 2
    table(data$MRSresid2cat)
    
    formula_text <- paste0("Surv(time = ", time_to, ", event = ", y, ") ~ MRSresid2cat")
    formula <- as.formula(formula_text)
    
    cox_model <- coxph(formula, data = data)
    
    coef_summary <- summary(cox_model)$coefficients
    
    results[i,] <- c(i,
                          coef_summary["MRSresid2cat", "exp(coef)"],
                          coef_summary["MRSresid2cat", "se(coef)"],
                          coef_summary["MRSresid2cat", "z"],
                          coef_summary["MRSresid2cat", "Pr(>|z|)"])
    
    ph_tests[[i]] <- cox.zph(cox_model)
    
  }  
  
  results_df <- results %>% as.data.frame()
  
  print(range(results_df$HR)) #[1] 0.6123488 1.5367513
  print(sum(results_df$P_value < 0.05)) #13)
  
  names(ph_tests) <- 1:1000
  violating_vars <- names(ph_tests)[sapply(ph_tests, function(x) x$table[1, 3] < 0.05)]
  print(violating_vars)
  
  return(results_df)
}

matrices <- as.list(mget(matrix_names))
names(matrices) <- dataset_names

chd_names <- list("CHD", "CHD", "chda")
names(chd_names) <- dataset_names
time_event_names <- list("days", "days", "chdatt")
names(time_event_names) <- dataset_names

cox_results <- lapply(dataset_names, function(dataset){
  cox_fun(
    resid_matrix = matrices[[dataset]],
    df = datasets[[dataset]],
    y = chd_names[[dataset]], 
    time_to = time_event_names[[dataset]]
  )
})
#V1
#0.6208869 1.5278398
#[1] 13
#[1] "224" "634" "678" "687" "850" "940"
#V2
# 0.5130779 1.2741019
# [1] 20
# NA
#MESA
# 0.6032669 2.6788420
# [1] 314
# [1] "472" "615" "785" "818" "824" "901"

names(cox_results) <- paste0(c("JHSv1", "JHSv2", "MESA"), "_cox_results")
list2env(cox_results, envir = .GlobalEnv)

# Save the list
saveRDS(cox_results, file = "cox_results_FDR.rds")

####### Permutation meta-analysis ########

JHSv1_cox_results$study <- ("EPICv1")
JHSv2_cox_results$study <- ("EPICv2")
MESA_cox_results$study <- ("MESA")

get_meta_results <- function(studies){
  
  v1 <- studies[[1]]
  v2 <- studies[[2]]
  mesa <- studies[[3]]
  
  # initialize results
  sitebind <- NULL
  run_meta <- NULL
  run_coef <- NULL
  meta_coef <- NULL
  i <- 1
  
  for (i in seq_along(v1$Iteration)){
    meta_cox <- rbind(v1[i,], v2[i,], mesa[i,])
    meta <- meta.summaries(log(HR), SE, method = c("fixed"), logscale = FALSE,
                           names = study, conf.level = 0.95, data = meta_cox,
                           subset = NULL)
    run_coef <- c(v1$Iteration[i], meta$summary, meta$se.summary, meta$het[3], meta$test, NA)
    meta_coef <- rbind(meta_coef, run_coef)
  }
  
  meta_coef_df <- as.data.frame(meta_coef)
  colnames(meta_coef_df) <- c("Iteration", "Coefficient", "SE", "p_het", "z", "p_meta", "HR")
  
  meta_coef_df$HR <- exp(meta_coef_df$Coefficient)
  
  meta_coef_df$q_meta <- p.adjust(meta_coef_df$p_meta, method = "fdr")
  meta_coef_df$bonferroni <- p.adjust(meta_coef_df$p_meta, method = "bonferroni")
  
  print(sum(meta_coef_df$p_het < 0.05)) # 98
  print(sum(meta_coef_df$p_meta < 0.05)) # 9
  print(sum(meta_coef_df$q_meta < 0.05)) # 0
  print(sum(meta_coef_df$bonferroni < 0.05)) # 0
  
  meta_coef_df_fixed <- meta_coef_df %>% filter(p_het > 0.05)
  meta_coef_df_fixed$method <- "FE"
  
  # Run random effects for heterogeneous permutations
  sig_hg <- meta_coef_df %>% filter(p_het < 0.05) %>% pull(Iteration)
  
  v1_hg <- v1[v1$Iteration %in% sig_hg,]
  v2_hg <- v2[v2$Iteration %in% sig_hg,]
  mesa_hg <- mesa[mesa$Iteration %in% sig_hg,]
  
  sitebind <- NULL
  run_meta <- NULL
  run_coef <- NULL
  meta_coef <- NULL
  i <- 1
  
  for (i in seq_along(v1_hg$Iteration)){
    meta_cox <- rbind(v1_hg[i,], v2_hg[i,], mesa_hg[i,])
    meta <- meta.summaries(log(HR), SE, method = c("random"), logscale = FALSE,
                           names = study, conf.level = 0.95, data = meta_cox,
                           subset = NULL)
    run_coef <- c(v1_hg$Iteration[i], meta$summary, meta$se.summary, meta$het[3], meta$test, NA)
    meta_coef <- rbind(meta_coef, run_coef)
  }
  
  meta_coef_df <- as.data.frame(meta_coef)
  colnames(meta_coef_df) <- c("Iteration", "Coefficient", "SE", "p_het", "z", "p_meta", "HR")
  
  meta_coef_df$HR <- exp(meta_coef_df$Coefficient)
  
  meta_coef_df$q_meta <- p.adjust(meta_coef_df$p_meta, method = "fdr")
  meta_coef_df$bonferroni <- p.adjust(meta_coef_df$p_meta, method = "bonferroni")
  
  print(sum(meta_coef_df$p_het < 0.05)) # 98
  print(sum(meta_coef_df$p_meta < 0.05)) # 0
  print(sum(meta_coef_df$q_meta < 0.05)) # 0
  print(sum(meta_coef_df$bonferroni < 0.05)) # 0
  
  meta_coef_df$method <- "RE"
  
  meta_coef_df_all <- bind_rows(meta_coef_df_fixed, meta_coef_df)
  
  return(meta_coef_df_all)
}

study_list <- list(JHSv1_cox_results, JHSv2_cox_results, MESA_cox_results)
meta_results <- get_meta_results(study_list)

write.csv(meta_results, file = "MRS_CHD_permutations_FDR_meta.csv")

########## Calculate p-value and plot Z-scores ###########
meta_results <- read.csv("MRS_CHD_permutations_FDR_meta.csv")

z_obs <- 1.708754218
z_perm <- c(meta_results$z)
p_perm <- (sum(z_perm >= z_obs) + 1) / (length(z_perm) + 1) #[1] 0.00109989
df <- data.frame(z = z_perm)

theme_set(theme_gray(base_family = "arial"))
update_geom_defaults("label", list(family = "arial"))
update_geom_defaults("text", list(family = "arial"))

z_plot <- ggplot(df, aes(x = z)) + 
  geom_density(fill = "lightgrey", alpha = 0.6, color = "darkgrey", linewidth = 1.2) +
  geom_vline(xintercept = z_obs, color = "#d73027", linetype = "dashed", linewidth = 1) +
  ggplot2::annotate("label", 
           x = z_obs,
           y = max(density(z_perm)$y),
           label = paste("Observed z =", round(z_obs, 2), "\n", "p-value =", signif(p_perm, digits = 4)),
           hjust = 0.5,
           vjust = 0.55,
           color = "black", 
           fill = "white",
           linewidth = 0,
           size = 6) +
  labs(
    title = "MRS<sub>841</sub> Permutation Z-score Distribution",
    x = "Z-score",
    y = "Density"
  ) + 
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linetype = "solid", linewidth = 0.5),
    plot.title = element_markdown(colour = "black", size = 22, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 19),
    axis.text.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 19),
    axis.text.y = element_text(colour = "black", size = 16),
    legend.title = element_text(colour = "black", size = 18),
    legend.text = element_text(colour = "black", size = 18),
    panel.grid.minor = element_blank(),
  )

z_plot

MRS841_z <- z_plot

ggsave("JHS_MESA_MRS_CHD_permutations_FDR.png", MRS841_z, width = 10, 
       height = 7, units = "in", dpi = 300)


########################################
############## MRS 841 #################
########################################

###### Residualize MRS #########

library(dplyr)
library(tidyverse)

null_scores_matrix1 <- readRDS("JHSv1_MRS_permutations_Bonferroni.rds")
null_scores_matrix2 <- readRDS("JHSv2_MRS_permutations_Bonferroni.rds")
null_scores_matrix_mesa <- readRDS("MESA_MRS_permutations_Bonferroni.rds")

dataset_names <- c("jhs_v1", "jhs_v2", "mesa")
datasets <- as.list(mget(dataset_names))
names(datasets) <- dataset_names

matrix_names <- c("null_scores_matrix1", "null_scores_matrix2",
                  "null_scores_matrix_mesa")
matrices <- as.list(mget(matrix_names))
names(matrices) <- dataset_names

subjid_names <- list("Sample_Name", "SUBJECT_ID", "TOEID")
names(subjid_names) <- dataset_names

null_matrices <- lapply(dataset_names, function(dataset){
  get_null_matrix(
    null_scores_matrix = matrices[[dataset]],
    sample_name = subjid_names[[dataset]],
    df = datasets[[dataset]]
  )
})

#[1] V1 dim before filtering: 1709 10000
#[1] TRUE
#[1] V1 dim after filtering: 1427 10000
#[1] V2 dim before filtering: 1687 10000
#[1] TRUE
#[1] V2 dim after filtering: 1323 10000
#[1] MESA dim before filtering: 870 10000
#[1] TRUE
#[1] MESA dim after filtering: 708 10000

names(null_matrices) <- matrix_names
list2env(null_matrices, envir = .GlobalEnv)

matrix_names <- c("null_scores_matrix1", "null_scores_matrix2",
                  "null_scores_matrix_mesa")
matrices <- as.list(mget(matrix_names))
names(matrices) <- dataset_names

cov_names <- c("jhs_cov", "jhs_cov",
               "mesa_cov")
cov <- mget(cov_names)
names(cov) <- dataset_names

residuals <- lapply(dataset_names, function(dataset){
  get_resid(
    null_scores_matrix = matrices[[dataset]],
    covariates = cov[[dataset]],
    df = datasets[[dataset]]
  )
})

names(residuals) <- paste0(dataset_names, "_resid")
list2env(residuals, envir = .GlobalEnv)

# Save the list
saveRDS(residuals, file = "perm_resid_data_Bonferroni.rds")


#### Determine association with CHD ####
library(survival)

# Load resid data
loaded_list <- readRDS("perm_resid_data_Bonferroni.rds")
dataset_names <- c("jhs_v1", "jhs_v2", "mesa")
names(loaded_list) <- paste0(dataset_names, "_resid")
list2env(loaded_list, envir = .GlobalEnv)

# JHS V1
jhs_v1_chd <- read.csv("/proj/azannas/projects/Hazel/JHS/Input/chd_data.csv", header = TRUE)
jhs_v2_chd <- read.csv("/proj/azannas/projects/Hazel/JHS/Input/chd_data_v2.csv", header = TRUE)
all_jhs_chd_dat <- read.csv("/proj/azannas/projects/JHS/incevtchd_2020_update_dbgap_id.csv")

jhs_v1_chd <- jhs_v1_chd %>%
  left_join(all_jhs_chd_dat, by = "SUBJECT_ID", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y"))

jhs_v2_chd <- jhs_v2_chd %>%
  left_join(all_jhs_chd_dat, by = "SUBJECT_ID", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y"))

dim(jhs_v1_chd) #n = 1554
dim(jhs_v2_chd) #n = 1499

mesa_chd <- filter(mesa, exall != 1)
mesa_chd <- mesa_chd[complete.cases(mesa_chd[ , mesa_cov]), ]
nrow(mesa_chd) # n = 707

jhs_v1_chd <- jhs_v1_chd[complete.cases(jhs_v1_chd[ , jhs_cov]), ]
nrow(jhs_v1_chd) # n = 1321

jhs_v2_chd <- jhs_v2_chd[complete.cases(jhs_v2_chd[ , jhs_cov]), ]
nrow(jhs_v2_chd) # n = 1227

dataset_names <- c("jhs_v1_chd", "jhs_v2_chd", "mesa_chd")
datasets <- as.list(mget(dataset_names))
names(datasets) <- dataset_names

matrix_names <- names(loaded_list)
matrices <- as.list(mget(matrix_names))
names(matrices) <- dataset_names

subjid_names <- list("Sample_Name", "SUBJECT_ID", "TOEID")
names(subjid_names) <- dataset_names

perm_resid_dat <- lapply(dataset_names, function(dataset){
  get_perm_resid_dat(
    resid_matrix = matrices[[dataset]],
    sample_name = subjid_names[[dataset]],
    df = datasets[[dataset]]
  )
})

#[1] TRUE
#[1] V1 dim after filtering: 1321 10000
#[1] TRUE
#[1] V2 dim after filtering: 1227 10000
#[1] TRUE
#[1] MESA dim after filtering: 707 10000

names(perm_resid_dat) <- matrix_names
list2env(perm_resid_dat, envir = .GlobalEnv)

matrices <- as.list(mget(matrix_names))
names(matrices) <- dataset_names

chd_names <- list("CHD", "CHD", "chda")
names(chd_names) <- dataset_names
time_event_names <- list("days", "days", "chdatt")
names(time_event_names) <- dataset_names

cox_results <- lapply(dataset_names, function(dataset){
  cox_fun(
    resid_matrix = matrices[[dataset]],
    df = datasets[[dataset]],
    y = chd_names[[dataset]], 
    time_to = time_event_names[[dataset]]
  )
})
#V1
# 0.4766675 1.8473034
# [1] 391
#"111"  "128"  "221"  "245"  "251"  "258"  "271"  "283"  "308"  "347"  "354"  "372" 
#[13] "386"  "394"  "449"  "457"  "477"  "502"  "521"  "528"  "599"  "608"  "641"  "650" 
#[25] "674"  "675"  "871"  "914"  "917"  "934"  "1000"
#V2
# 0.3391732 2.0628210
# [1] 481
# "19"  "20"  "41"  "71"  "139" "157" "213" "230" "433" "470" "526" "557" "581" "586"
# [15] "623" "795" "822" "834" "846" "859" "901" "907" "985" "990"
#MESA
# 0.290107 4.453309
# [1] 949
# [1] "10"  "18"  "20"  "21"  "28"  "41"  "60"  "91"  "107" "113" "122" "127" "131" "145"
# [15] "151" "156" "159" "162" "178" "183" "216" "225" "228" "232" "237" "250" "272" "276"
# [29] "283" "302" "318" "320" "334" "341" "352" "356" "381" "418" "438" "482" "496" "501"
# [43] "522" "534" "537" "582" "586" "589" "609" "622" "626" "636" "663" "666" "669" "713"
# [57] "724" "746" "766" "781" "793" "805" "814" "817" "819" "830" "848" "851" "852" "854"
# [71] "892" "921" "940" "941" "957" "961" "977" "994" "998"

names(cox_results) <- paste0(c("JHSv1", "JHSv2", "MESA"), "_cox_results")
list2env(cox_results, envir = .GlobalEnv)

# Save the list
saveRDS(cox_results, file = "cox_results_Bonferroni.rds")

####### Permutation meta-analysis #########

library(rmeta)
library(tidyverse)
library(dplyr)

JHSv1_cox_results$study <- ("EPICv1")
JHSv2_cox_results$study <- ("EPICv2")
MESA_cox_results$study <- ("MESA")

study_list <- list(JHSv1_cox_results, JHSv2_cox_results, MESA_cox_results)
meta_results <- get_meta_results(study_list)
# [1] 679
# [1] 510
# [1] 0
# [1] 0
# [1] 679
# [1] 1
# [1] 0
# [1] 0

write.csv(meta_results, file = "MRS_CHD_permutations_Bonferroni_meta.csv")

######## Calculate p-value and plot Z-scores ##########
meta_results <- read.csv("MRS_CHD_permutations_Bonferroni_meta.csv")

z_obs <- 2.09374772
z_perm <- c(meta_results$z)
p_perm <- (sum(z_perm >= z_obs) + 1) / (length(z_perm) + 1) #[1] 0.00869913

df <- data.frame(z = z_perm)

z_plot <- ggplot(df, aes(x = z)) + 
  geom_density(fill = "lightgrey", alpha = 0.6, color = "darkgrey", linewidth = 1.2) +
  geom_vline(xintercept = z_obs, color = "#d73027", linetype = "dashed", linewidth = 1) +
  ggplot2::annotate("label", 
                    x = z_obs,
                    y = max(density(z_perm)$y),
                    label = paste("Observed z =", round(z_obs, 2), "\n", "p-value =", signif(p_perm, digits = 2)),
                    hjust = 0.5,
                    vjust = 0.55,
                    color = "black", 
                    fill = "white",
                    linewidth = 0,
                    size = 6) +
  labs(
    title = "MRS<sub>13</sub> Permutation Z-score Distribution",
    x = "Z-score",
    y = "Density"
  ) + 
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),,
    panel.grid.major = element_line(colour = "grey90", linetype = "solid", linewidth = 0.5),
    plot.title = element_markdown(colour = "black", size = 22, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(colour = "black", size = 19),
    axis.text.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 19),
    axis.text.y = element_text(colour = "black", size = 16),
    legend.title = element_text(colour = "black", size = 18),
    legend.text = element_text(colour = "black", size = 18),
    panel.grid.minor = element_blank()
  )

z_plot
MRS13_z <- z_plot

ggsave("JHS_MESA_MRS_CHD_permutations_Bonferroni.png", MRS13_z, width = 10, 
       height = 7, units = "in", dpi = 300)

MRS841_z <- MRS841_z + theme(plot.tag = element_text(size = 22, face = "bold", family = "arial"))
MRS13_z  <- MRS13_z  + theme(plot.tag = element_text(size = 22, face = "bold", family = "arial"))

combined_plot <- MRS841_z + MRS13_z +
  plot_annotation(tag_levels = list(c("C.", "D.")))

combined_plot

ggsave("JHS_MESA_MRS_CHD_permutations_combined_plot.png", combined_plot, width = 20, 
       height = 7, units = "in", dpi = 300)
