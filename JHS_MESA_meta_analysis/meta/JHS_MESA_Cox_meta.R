####################################################
########## COX REGRESSION META-ANALYSIS ###########
####################################################

library(rmeta)

####### FDR-adjusted sites ######
mesa_path <- "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/002 Zannas_lab/MRS TOPMed/MRS_WHI_MESA/Output/" # directory for MESA Cox results
mesa <- read.csv(paste0(mesa_path, "MESA_cox_results_fdr.csv"), header = TRUE)

jhs_path <- "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/002 Zannas_lab/MRS TOPMed/JHS_MESA_meta_analysis/Input/"
v1 <- read.csv(paste0(jhs_path, "v1_cox_results_fdr.csv"), header = TRUE)
v2 <- read.csv(paste0(jhs_path, "v2_cox_results_fdr.csv"), header = TRUE)

v1$study <- ("JHSv1")
v2$study <- ("JHSv2")
mesa$study <- ("MESA")
meta_cox <- rbind(v1, v2, mesa)
meta_cox <- as.data.frame(meta_cox)
glimpse(meta_cox)

meta <- meta.summaries(log(HR), SE, method = c("fixed"), logscale = FALSE,
                       names = study, conf.level = 0.95, data = meta_cox,
                       subset = NULL)
run_coef <- c(meta$summary, meta$se.summary, meta$het[3], meta$test, NA)
names(run_coef) <- c("Coefficient", "SE", "p_het", "z", "p_meta", "HR")
print(run_coef)
run_coef[6] <- exp(run_coef[1])
print(meta) # 95% CI (-0.0346, 0.505)

# Compute CI: b - (1.96 * SE), b + (1.96*SE)
run_coef$SE

computeCI <- function(hr, se){
  
  HR <- hr
  SE <- se
  
  lowerCI <- HR - (1.96 * SE)
  upperCI <- HR + (1.96 * SE)
  
  print(paste0("(", round(lowerCI, 2), ", ", round(upperCI, 2), ")"))
}

ci_vals <- computeCI(run_coef["HR"], run_coef["SE"])

run_coef <- run_coef %>% t() %>% as.data.frame()

rownames(run_coef) <- "CHD_FDR"
run_coef$CI <- ci_vals

CHD_FDR <- run_coef


meta_cox_fdr_all <- run_coef
setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Output")
write.csv(meta_cox_fdr_all, "JHS_MESA_CHD_MRS_FDR.csv")

####### Bonferroni-adjusted sites ######
mesa <- read.csv(paste0(mesa_path, "MESA_cox_results_bonf.csv"), header = TRUE)

v1 <- read.csv(paste0(jhs_path, "v1_cox_results_bonf.csv"), header = TRUE)
v2 <- read.csv(paste0(jhs_path, "v2_cox_results_bonf.csv"), header = TRUE)

v1$study <- ("JHSv1")
v2$study <- ("JHSv2")
mesa$study <- ("MESA")
meta_cox <- rbind(v1, v2, mesa)
meta_cox <- as.data.frame(meta_cox)
glimpse(meta_cox)

meta <- meta.summaries(log(HR), SE, method = c("fixed"), logscale = FALSE,
                       names = study, conf.level = 0.95, data = meta_cox,
                       subset = NULL)
run_coef <- c(meta$summary, meta$se.summary, meta$het[3], meta$test, NA)
names(run_coef) <- c("Coefficient", "SE", "p_het", "z", "p_meta", "HR")
print(run_coef)
run_coef[6] <- exp(run_coef[1])
print(meta) # 95% CI (0.0185, 0.561)

ci_vals <- computeCI(run_coef["HR"], run_coef["SE"])

run_coef <- run_coef %>% t() %>% as.data.frame()

rownames(run_coef) <- "CHD_Bonf"
run_coef$CI <- ci_vals

CHD_Bonf <- run_coef

meta_cox_bonf_all <- run_coef
setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis/Output")
write.csv(meta_cox_bonf_all, "JHS_MESA_CHD_MRS_Bonferroni.csv")

###############################################
########## PLOT KAPLAN-MEIER CURVES ###########
###############################################

############ FDR #############

mesa_fdr_all # See "MESA_Cox.R" for data prep
v1_fdr_all # See "JHS_Cox.R" for data prep
v2_fdr_all # See "JHS_Cox.R" for data prep

v1 <- subset(v1_fdr_all, select = c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex"))
v2 <- subset(v2_fdr_all, select = c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex"))
mesa <- subset(mesa_fdr_all, select = c("sidno", "chda", "chdatt", "MRSresid2cat", "gender1"))

colnames(mesa) <- c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex")
data <- rbind(v1, v2, mesa)

surv_obj <- Surv(data$days_to_event, data$CHD == 1)
km_fit <- survfit(surv_obj ~ MRSresid2cat, data = data)

xmax <- range(data$days_to_event, na.rm = TRUE)[2]

library(ggtext)
event_plot <- ggsurvplot(km_fit, 
                         data = data, 
                         pval = FALSE, 
                         conf.int = TRUE,
                         fun = "event",  # This inverts to show event probability
                         legend.labs = c("Lower", "Higher"),
                         legend.title = "MRS Score",
                         title = "CHD Cumulative Incidence by MRS<sub>841</sub>",
                         xlab = "Time to Event (Days)",
                         ylab = "CHD Probability",
                         xlim = c(0, xmax))$plot + 
  theme(panel.grid.major = element_line(colour = "grey90", linetype = "solid", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linetype = "solid", linewidth = 0.5),
        plot.title = element_markdown(colour = "black", size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 19),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 19),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  annotate("text", 
         x = 2000, # x-coordinate for the text
         y = 0.1, # y-coordinate for the text
         label = "p-val = 0.0875",
         size = 6,
         color = "black")

event_plot
setwd("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/002 Zannas_lab/MRS TOPMed/JHS_MESA_meta_analysis/Results/")
ggsave("MRS_CHD_FDR.png", event_plot, scale = 1, width = 10, height = 7, units = "in")

############ Bonferroni #############

# Prepare data using lines 142 - 187 in "Cox_MESA_MRS_resid.R"
mesa_bonf_all # See "MESA_Cox.R" for data prep
v1_bonf_all # See "JHS_Cox.R" for data prep
v2_bonf_all # See "JHS_Cox.R" for data prep

#### All ####
v1 <- subset(v1_bonf_all, select = c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex"))
v2 <- subset(v2_bonf_all, select = c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex"))
mesa <- subset(mesa_bonf_all, select = c("sidno", "chda", "chdatt", "MRSresid2cat", "gender1"))

colnames(mesa) <- c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex")
data <- rbind(v1, v2, mesa)

surv_obj <- Surv(data$days_to_event, data$CHD == 1)
km_fit <- survfit(surv_obj ~ MRSresid2cat, data = data)

event_plot <- ggsurvplot(km_fit, 
                         data = data, 
                         pval = FALSE, 
                         conf.int = TRUE,
                         fun = "event",  # This inverts to show event probability
                         legend.labs = c("Lower", "Higher"),
                         legend.title = "MRS Score",
                         title = "CHD Cumulative Incidence by MRS<sub>13</sub>",
                         xlab = "Time to Event (Days)",
                         ylab = "CHD Probability",
                         xlim = c(0, xmax))$plot + 
  theme(panel.grid.major = element_line(colour = "grey90", linetype = "solid", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linetype = "solid", linewidth = 0.5),
        plot.title = element_markdown(colour = "black", size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 19),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 19),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  annotate("text", 
           x = 2000,
           y = 0.1,
           label = "p-val = 0.0363",
           size = 6,
           color = "black")

event_plot

setwd("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/002 Zannas_lab/MRS TOPMed/JHS_MESA_meta_analysis/Results/")
ggsave("MRS_CHD_Bonferroni.png", event_plot, width = 10, height = 7, units = "in")
