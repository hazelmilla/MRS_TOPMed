library(survival)
library(survminer)
library(tidyverse)
library(dplyr)

dir = "/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/V1/"
setwd(dir)

data_unfiltered <- read.csv("Input_files/chd_data.csv", header = TRUE)
dim(data_unfiltered)
#[1] 1554   40

####### FDR THRESHOLD ########

# Filter data
data <- data_unfiltered %>% drop_na(CHD) %>% drop_na(days_to_event) %>%
  drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>% drop_na(MRS2cat)
dim(data)
#[1] 1291   40

# Ensure MRS with median split is a factor with proper labels
data$MRS2cat <- factor(data$MRS2cat, 
                             levels = c(1, 2), 
                             labels = c("Low MRS", "High MRS"))

##### Fit Cox proportional hazards model #####

cox_model <- coxph(Surv(days_to_event, CHD) ~ MRS2cat + sex + age + 
                     currentSmoker + BMI + Income + edu3cat + alc + 
                     marital2cat + nbSESanascore + 
                     PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                     PC8 + PC9 + PC10 + NK + Mono + Gran + 
                     Bcell + CD8T + CD4T, 
                   data = data)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "OR", "SE",
                            "z", "P_val")

setwd(outputdir)
write.csv(coefficients, file = "Output_files/v1_cox_results_fdr.csv")

# check assumptions
test_ph <- cox.zph(cox_model)
print(test_ph)
# if p-value is significant, stratify by that variable
# no p values are significant

ggforest(cox_model, data = data) # visualize results

####### Generate Kaplan-Meier curve #########

setwd(inputdir)
# Create survival object
surv_obj <- Surv(data$days_to_event, data$CHD == 1)

# Fit Kaplan-Meier survival curves
km_fit <- survfit(surv_obj ~ MRS2cat, data = data)

survival_plot <- ggsurvplot(km_fit, 
                            data = data, 
                            pval = TRUE, 
                            conf.int = TRUE,
                            legend.labs = c("Low", "High"),
                            legend.title = "MRS Level",
                            title = "CHD Risk by MRS level",
                            xlab = "Time to Event in Days",
                            ylab = "CHD Probability")$plot

survival_plot

event_plot <- ggsurvplot(km_fit, 
                         data = data, 
                         pval = TRUE, 
                         conf.int = TRUE,
                         fun = "event",  # This inverts to show event probability
                         legend.labs = c("Low", "High"),
                         legend.title = "MRS Level",
                         title = "CHD Cumulative Incidence by MRS level",
                         xlab = "Time to Event in Days",
                         ylab = "CHD Probability")$plot

event_plot

###### BONFERRONI THRESHOLD #########
data <- data_unfiltered %>% drop_na(CHD) %>% drop_na(days_to_event) %>%
  drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>% drop_na(MRS2cat_bonferroni)
dim(data)
#[1] 1291   39

# Ensure MRS_med_split is a factor with proper labels
data$MRS2cat_bonferroni <- factor(data$MRS2cat_bonferroni, 
                           levels = c(1, 2), 
                           labels = c("Low MRS", "High MRS"))

##### Fit Cox proportional hazards model #####
cox_model <- coxph(Surv(days_to_event, CHD) ~ MRS2cat_bonferroni + sex + age + 
                     currentSmoker + BMI + Income + edu3cat + alc + 
                     marital2cat + nbSESanascore +
                     PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                     PC8 + PC9 + PC10 + NK + Mono + Gran + 
                     Bcell + CD8T + CD4T, 
                   data = data)

# check assumptions of the model
test_ph <- cox.zph(cox_model)
print(test_ph)

print(ifelse(test_ph$table[,3] < 0.05, "yes", "no")[test_ph$table[,3] < 0.05])
# if p-value is significant, you need to stratify by that variable
# no p values are significant

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "OR", "SE",
                            "z", "P_val")

write.csv(coefficients, file = "Output_files/v1_cox_results_bonferroni.csv")

ggforest(cox_model, data = data) # visualize results

####### Generate Kaplan-Meier curve #########

# Create survival object
surv_obj <- Surv(data$days_to_event, data$CHD == 1)

# Fit Kaplan-Meier survival curves
km_fit <- survfit(surv_obj ~ MRS2cat_bonferroni, data = data)

survival_plot <- ggsurvplot(km_fit, 
                            data = data, 
                            pval = TRUE, 
                            conf.int = TRUE,
                            legend.labs = c("Low", "High"),
                            legend.title = "MRS Level",
                            title = "CHD Risk by MRS level",
                            xlab = "Time to Event in Days",
                            ylab = "CHD Probability")$plot

survival_plot

event_plot <- ggsurvplot(km_fit, 
                         data = data, 
                         pval = TRUE, 
                         conf.int = TRUE,
                         fun = "event",  # This inverts to show event probability
                         legend.labs = c("Low", "High"),
                         legend.title = "MRS Level",
                         title = "CHD Cumulative Incidence by MRS level",
                         xlab = "Time to Event in Days",
                         ylab = "CHD Probability")$plot

event_plot
