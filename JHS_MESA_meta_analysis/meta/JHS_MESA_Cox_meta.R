#


###############################################
########## PLOT KAPLAN-MEIER CURVES ###########
###############################################

######## FDR #############
# Prep MESA data using lines 5 - 52 "Cox_MESA_MRS_resid.R"
mesa_fdr_all
mesa_fdr_f #71-89
mesa_fdr_m #106-124

# Prep JHS v1 data using lines 6 - 53 "MRSresidual_CHD_JHS.R"
v1_fdr_all
v1_fdr_f #73-99
v1_fdr_m #117-143

# Prep JHS v2 data using lines 164 - 209 "MRSresidual_CHD_JHS.R"
v2_fdr_all
v2_fdr_f #230-256
v2_fdr_m #274-300

#### All ####
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

#### Female ####
v1 <- subset(v1_fdr_f, select = c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex"))
v2 <- subset(v2_fdr_f, select = c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex"))
mesa <- subset(mesa_fdr_f, select = c("sidno", "chda", "chdatt", "MRSresid2cat", "gender1"))

colnames(mesa) <- c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex")
data_f <- rbind(v1, v2, mesa)

surv_obj <- Surv(data_f$days_to_event, data_f$CHD == 1)
km_fit <- survfit(surv_obj ~ MRSresid2cat, data = data_f)

event_plot <- ggsurvplot(km_fit, 
                         data = data_f, 
                         pval = FALSE, 
                         conf.int = TRUE,
                         fun = "event",  # This inverts to show event probability
                         legend.labs = c("Low", "High"),
                         legend.title = "MRS Level",
                         title = "CHD Cumulative Incidence by MRS Level (FDR) - Female",
                         xlab = "Time to Event in Days",
                         ylab = "CHD Probability")$plot + 
  theme(panel.grid.major = element_line(colour = "grey90", linetype = "solid", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linetype = "solid", linewidth = 0.5)) +
  annotate("text", 
           x = 2000,
           y = 0.1,
           label = "p = 0.4047",
           size = 5,
           color = "black")

event_plot
ggsave("MRS_CHD_FDR_female.png", event_plot, width = 10, height = 7, units = "in")

#### Male ####
v1 <- subset(v1_fdr_m, select = c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex"))
v2 <- subset(v2_fdr_m, select = c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex"))
mesa <- subset(mesa_fdr_m, select = c("sidno", "chda", "chdatt", "MRSresid2cat", "gender1"))

colnames(mesa) <- c("SUBJECT_ID", "CHD", "days_to_event", "MRSresid2cat", "sex")
data_m <- rbind(v1, v2, mesa)

surv_obj <- Surv(data_m$days_to_event, data_m$CHD == 1)
km_fit <- survfit(surv_obj ~ MRSresid2cat, data = data_m)

event_plot <- ggsurvplot(km_fit, 
                         data = data_m, 
                         pval = FALSE, 
                         conf.int = TRUE,
                         fun = "event",  # This inverts to show event probability
                         legend.labs = c("Low", "High"),
                         legend.title = "MRS Level",
                         title = "CHD Cumulative Incidence by MRS Level (FDR) - Male",
                         xlab = "Time to Event in Days",
                         ylab = "CHD Probability")$plot + 
  theme(panel.grid.major = element_line(colour = "grey90", linetype = "solid", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linetype = "solid", linewidth = 0.5)
  ) +
  annotate("text", 
           x = 2000,
           y = 0.1,
           label = "p = 0.5781",
           size = 5,
           color = "black")

event_plot

ggsave("MRS_CHD_FDR_male.png", event_plot, width = 10, height = 7, units = "in")


######## Bonferroni #############

# Prepare data using lines 142 - 187 in "Cox_MESA_MRS_resid.R"
mesa_bonf_all
mesa_bonf_f #206-223
mesa_bonf_m #240-257

# Prepare data using lines 326 - 372 in "MRSresidual_CHD_JHS.R"
v1_bonf_all
v1_bonf_f #391-417
v1_bonf_m #435-460

# Prepare data using lines 481 - 525 in "MRSresidual_CHD_JHS.R"
v2_bonf_all
v2_bonf_f # 546-571
v2_bonf_m #589-614

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
