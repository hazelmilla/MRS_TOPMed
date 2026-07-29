# Run Cox proportional hazards regression for all CHD events

library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggtext)
library(survminer)
library(survival)
library(rmeta)
library(patchwork)

############################
####### Prepare data #######
############################

dir = "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/000 Zannas_lab/MRS TOPMed/JHS_MESA_meta_analysis/Cox_update_June2026/"
setwd(dir)

jhs_v1 <- read.csv("chd_data_v1.csv", header = TRUE)
jhs_v2 <- read.csv("chd_data_v2.csv", header = TRUE)
mesa <- read.csv("prepared_data_MESA.csv", header = TRUE)

# JHS
jhs_pheno <- read.csv("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/JHS/Pheno_Data/analysis1_2020_update_dbgap_id.csv",
                      header = TRUE)
jhs_pheno <- subset(jhs_pheno, select = c("SUBJECT_ID", "everSmoker"))

jhs_v1 <- left_join(jhs_v1, jhs_pheno, by = "SUBJECT_ID")
jhs_v2 <- left_join(jhs_v2, jhs_pheno, by = "SUBJECT_ID")
nrow(jhs_v1) # n = 1554
nrow(jhs_v2) # n = 1499

# MESA 
# Cigarette smoking status - recode cig1c -> lifesmk1 0 as 0, 1 and 2 as 1

# Get all pheno data 
mesa_pheno <- read.csv("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/MESA/SHAReE1_ZannasPheno_20230623.csv",
         header = TRUE)
mesa_pheno <- subset(mesa_pheno, select = c("sidno", "cig1c"))
mesa_pheno <- mutate(mesa_pheno, lifesmk1 = ifelse(cig1c == 2, 1, cig1c))

mesa <- left_join(mesa, mesa_pheno, by = "sidno")
nrow(mesa) # n = 870


# Get all CHD info
jhs_events <- read.csv("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/JHS/Pheno_Data/incevtchd_2020_update_dbgap_id.csv",
                       header = TRUE)

jhs_v1 <- jhs_v1 %>%
  left_join(jhs_events, by = "SUBJECT_ID", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y"))

jhs_v2 <- jhs_v2 %>%
  left_join(jhs_events, by = "SUBJECT_ID", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y"))
            
dim(jhs_v1) #n = 1554
dim(jhs_v2) #n = 1499

# get hard and MI event data
mesa_hard_mi <- read.csv("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/MESA/SHAReE1_ZannasEV_20230623.csv",
                         header = TRUE)

mesa <- mesa %>%
  left_join(mesa_hard_mi, by = "sidno", suffix = c("", ".y")) %>%
  dplyr::select(-ends_with(".y"))
dim(mesa) #n = 870

mesa <- filter(mesa, exall != 1) # Filter out prior history of CHD from MESA
nrow(mesa) # n = 869
# Prior history of CHD has already been removed for JHS

##################################################
####### Run Cox for each individual cohort #######
##################################################

# Define functions

remove_cov_na <- function(df, covariates){
  data <- as.data.frame(df)[complete.cases(as.data.frame(df)[ , c("MRS", covariates)]), ]
  return(data)
}

# Generate resid function
get_resid <- function(x, covariates, df){
  
  data <- as.data.frame(df)
  
  formula_text <- paste0(x, "~",
                         paste(c(covariates), collapse = " + "))
  formula <- as.formula(formula_text)
  residual <- lm(formula, data = data)
  
  data$residual <- resid(residual)
  
  med <- summary(data$residual)[3]
  data$Resid2cat[data$residual <= med] <- 1
  data$Resid2cat[data$residual > med] <- 2
  table(data$Resid2cat)
  
  resid_name <- paste0(x, "_resid2cat")
  
  idx <- which(colnames(data) == "Resid2cat")
  
  colnames(data)[idx] <- resid_name
  
  data$residual = NULL
  
  return(data)
}

cox_fun <- function(x, y, time_to, df){
  
  data <- as.data.frame(df)[complete.cases(as.data.frame(df)[ , c(x, y, time_to)]), ]
  print(nrow(data)) # Get n-values
  
  formula_text <- paste0("Surv(time = ", time_to, ", event = ", y, ") ~", x)
  formula <- as.formula(formula_text)
  print(formula)
  
  model <- coxph(formula, data = data)

  print(cox.zph(model)) # see assumption violation p-values

  # Summarize results
  coefficients <- as.data.frame(summary(model)$coefficients)

  # get 95% confidence interval
  ci <- summary(model)$conf.int

  # get coefficients
  coefficients$lower95ci <- ci[3]
  coefficients$upper95ci <- ci[4]
  colnames(coefficients) <- c("Estimate", "HR", "SE",
                              "z", "P_val", "Lower CI", "Upper CI")
  return(coefficients)
}

# Define CI function 
# Compute CI: b - (1.96 * SE), b + (1.96*SE)
computeCI <- function(hr, se){
  
  HR <- hr
  SE <- se
  
  lowerCI <- HR - (1.96 * SE)
  upperCI <- HR + (1.96 * SE)
  
  print(paste0("(", signif(lowerCI, 5), ", ", signif(upperCI, digits = 5), ")"))
}

# Define meta-analysis function
get_meta_result <- function(df){
  
  meta <- meta.summaries(log(HR), SE, method = c("fixed"), logscale = FALSE,
                         names = study, conf.level = 0.95, data = df,
                         subset = NULL)
  run_coef <- c(meta$summary, meta$se.summary, meta$het[3], meta$test, NA)
  names(run_coef) <- c("Coefficient", "SE", "p_het", "z", "p_meta", "HR")
  run_coef[6] <- exp(run_coef[1])
  print(run_coef)
  
  ci_vals <- computeCI(run_coef["HR"], run_coef["SE"])
  
  run_coef <- run_coef %>% t() %>% as.data.frame()
  run_coef$CI <- ci_vals
  
  return(run_coef)
}

# Define covariates
# JHS
jhs_cov <- c("sex", "age", "currentSmoker", "BMI", "Income",
             "edu3cat", "alc", "marital2cat", "nbSESanascore",
             "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",
             "PC9", "PC10", "NK", "Mono", "Gran", "Bcell", "CD8T",
             "CD4T")
# MESA 
mesa_cov <- c("gender1", "AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
              "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
              "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", 
              "Sentrix_ID", "Sample_Well", "Sample_Plate", "NK", "Gran", "Mono", 
              "Bcell", "CD8T", "CD4T")


# Remove NA covariates from JHS and MESA


# prepare data for remove NA function
cov_list <- list(jhs_cov, jhs_cov, mesa_cov)
datasets = list(jhs_v1, jhs_v2, mesa)
names(datasets) = c("jhs_v1", "jhs_v2", "mesa")
names(cov_list) <- names(datasets)

data_no_na <- lapply(names(datasets), function(set){
  remove_cov_na(df = datasets[[set]], 
  covariates = c(unlist(cov_list[[set]])))
})

names(data_no_na) <- paste0(names(datasets), "_chd") # Rename datasets
list2env(data_no_na, envir = .GlobalEnv) # move from list to environment

# Get n values for all CHD events
sum(!is.na(jhs_v1_chd$CHD)) # current smoking n = 1291; ever smoking n = 1299
sum(!is.na(jhs_v2_chd$CHD)) # current smoking n = 1202, ever smoking n = 1211
sum(!is.na(mesa_chd$chda)) # current smoking n = 707; ever smoking n = 707

# Get MRS dichotomized residuals for MRS_841 (FDR)
# prepare data for get_resid function
dataset_names <- names(data_no_na)
datasets <- mget(dataset_names)
names(datasets) <- dataset_names
names(cov_list) <- names(datasets)

fdr_residuals <- lapply(names(datasets), function(set){
  get_resid(x = "MRS", 
            covariates = c(unlist(cov_list[[set]])), 
            df = datasets[[set]])
})

names(fdr_residuals) <- names(datasets) 
list2env(fdr_residuals, envir = .GlobalEnv) # move updated datasets to environment

# Get MRS dichotomized residuals for MRS_13 (Bonferroni)
# prepare data for get_resid function
colnames(mesa_chd)[which(colnames(mesa_chd) == "MRS_bonf")] <- "MRS_bonferroni"
datasets <- mget(dataset_names)

bonf_residuals <- lapply(names(datasets), function(set){
  get_resid(x = "MRS_bonferroni", 
            covariates = c(unlist(cov_list[[set]])), 
            df = datasets[[set]])
})

names(bonf_residuals) <- names(datasets)
list2env(bonf_residuals, envir = .GlobalEnv)

# Run Cox regression
datasets <- mget(dataset_names)
names(datasets) <- dataset_names

chd_names <- list("CHD", "CHD", "chda")
names(chd_names) <- names(datasets)
time_event_names <- list("days", "days", "chdatt")
names(time_event_names) <- names(datasets)

cox_results <- lapply(names(datasets), function(set){
  cox_fun(
    x = "MRS_resid2cat", 
    y = chd_names[[set]], 
    time_to = time_event_names[[set]], 
    df = datasets[[set]])
})

fdr_cox_results <- bind_rows(cox_results)
rownames(fdr_cox_results) <- names(datasets)

cox_results <- lapply(names(datasets), function(set){
  cox_fun(
    x = "MRS_bonferroni_resid2cat", 
    y = chd_names[[set]], 
    time_to = time_event_names[[set]], 
    df = datasets[[set]])
})

bonf_cox_results <- bind_rows(cox_results)
rownames(bonf_cox_results) <- names(datasets)

# Current smoking Cox assumption violation p-vals: 
# FDR - v1 p = 0.43, v2 = 0.45, mesa = 0.33
# Bonferroni - v1 = 0.48, v2 = 0.032, mesa = 0.59

############################################
####### Meta-analyze across cohorts ########
############################################

# FDR
fdr_cox_results$study <- c("JHSv1", "JHSv2", "MESA")
glimpse(fdr_cox_results)

bonf_cox_results$study <- c("JHSv1", "JHSv2", "MESA")
glimpse(bonf_cox_results)

results <- list(fdr_cox_results, bonf_cox_results)
names(results) <- c("fdr_cox_results", "bonf_cox_results")
meta_results <- lapply(names(results), function(res){
    get_meta_result(df = results[[res]])
    })

meta_results <- bind_rows(meta_results)
rownames(meta_results) = names(results)

write.csv(meta_results, "JHS_MESA_all_CHD_meta.csv")

##################################
####### Generate km curves ########
##################################

theme_set(theme_gray(base_family = "arial"))
update_geom_defaults("label", list(family = "arial"))
update_geom_defaults("text", list(family = "arial"))

v1 <- subset(jhs_v1_chd, select = c("SUBJECT_ID", "CHD", "days", "MRS_resid2cat", "MRS_bonferroni_resid2cat"))
v2 <- subset(jhs_v2_chd, select = c("SUBJECT_ID", "CHD", "days", "MRS_resid2cat", "MRS_bonferroni_resid2cat"))
mesa1 <- subset(mesa_chd, select = c("sidno", "chda", "chdatt", "MRS_resid2cat", "MRS_bonferroni_resid2cat"))

colnames(mesa1) <- c("SUBJECT_ID", "CHD", "days", "MRS_resid2cat", "MRS_bonferroni_resid2cat")
data <- rbind(v1, v2, mesa1)

data <- mutate(data, years_to_event = (days / 365)) # convert to years
xmax <- range(data$years_to_event, na.rm = TRUE)[2] # set max years for plot config

plot_configs <- list(
  list(group_var = "MRS_resid2cat",
       subscript = "841",
       pval_label = paste0("p-val = ", signif(meta_results$p_meta[1], digits = 2)),
       outfile   = "MRS_CHD_FDR.png"),
  
  list(group_var = "MRS_bonferroni_resid2cat",
       subscript = "13",
       pval_label = paste0("p-val = ", signif(meta_results$p_meta[2], digits = 2)),
       outfile   = "MRS_CHD_Bonf.png")) # list plot settings

make_km_plot <- function(cfg, data, xmax) {
  
  surv_formula <- as.formula(
    paste0("Surv(years_to_event, CHD == 1) ~ ", cfg$group_var)
  )
  
  # Force substitution of the actual formula object into the call,
  km_fit <- eval(bquote(survfit(.(surv_formula), data = data)))
  
  title_text <- paste0("CHD Cumulative Incidence by MRS<sub>", cfg$subscript, "</sub>")
  
  max_event_prob <- max(1 - km_fit$surv, na.rm = TRUE)
  y_center <- max_event_prob / 2
  
  event_plot <- ggsurvplot(km_fit,
                           data = data,
                           pval = FALSE,
                           conf.int = TRUE,
                           fun = "event",
                           legend.labs = c("Lower", "Higher"),
                           legend.title = "MRS Score",
                           title = title_text,
                           xlab = "Time to Event (Years)",
                           ylab = "CHD Probability",
                           xlim = c(0, xmax))$plot +
    theme(panel.grid.major = element_line(colour = "grey90", linetype = "solid", linewidth = 0.5),
          panel.grid.minor = element_line(color = "grey90", linetype = "solid", linewidth = 0.5),
          plot.title = element_markdown(colour = "black", size = 22, hjust = 0.5, face = "bold"),
          axis.title.x = element_text(size = 19),
          axis.text.x = element_text(size = 16),
          axis.title.y = element_text(size = 19),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18)) +
    annotate("text",
             x = 5.5,
             y = 1.6*y_center,
             label = cfg$pval_label,
             size = 6,
             color = "black")
  
  print(event_plot)
  ggsave(cfg$outfile, event_plot, scale = 1, width = 10, height = 7, units = "in")
  
  return(event_plot)
}

chd_plots <- lapply(plot_configs, make_km_plot, data = data, xmax = xmax) # generate plots

combined_plot <- chd_plots[[1]] + chd_plots[[2]] +
  plot_annotation(tag_levels = list(c("A.", "B."))) & 
  theme(plot.tag = element_text(size = 22, face = "bold")) # combine plots with patchwork

combined_plot # check combined plot

# save combined plot
ggsave("chd_km_curves.PNG", combined_plot, scale = 1, width = 20, height = 7, units = "in")

