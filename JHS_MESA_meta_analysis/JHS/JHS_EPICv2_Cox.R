# Cox regression - v2 data
library(survival)
library(survminer)
library(tidyverse)
library(dplyr)
library(ggplot2)

path = "/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/V2"
setwd(path)

# Fit a Cox proportional hazards model
# Plot survival time based on MRS & covariates

#######################################
###### Threshold 1 (FDR) #######
#######################################

data_unfiltered <- read.csv("chd_data_v2.csv", header = TRUE)

data <- data_unfiltered %>% drop_na(MRS2cat) %>% drop_na(CHD) %>% drop_na(days_to_event) %>%
  drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore)
dim(data)
#[1] 1202   43

# Ensure MRS_med_split is a factor with proper labels
data$MRS2cat <- factor(data$MRS2cat, 
                       levels = c(1, 2), 
                       labels = c("Low MRS", "High MRS"))

data$BMI3cat <- ntile(data$BMI, 3) # Split into tertiles due to assumption violation (p < 0.05)
data$NK3cat <- ntile(data$NK, 3) # Split into tertiles due to assumption violation (p < 0.05)
data$Gran3cat <- ntile(data$Gran, 3) # Split into tertiles due to assumption violation (p < 0.05)

######## Fit Cox proportional hazards regression model #########
cox_model <- coxph(Surv(days_to_event, CHD) ~ MRS2cat + sex + age + 
                     currentSmoker + strata(BMI3cat) + Income + edu3cat + alc + 
                     marital2cat + nbSESanascore + 
                     PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                     PC8 + PC9 + PC10 + strata(NK3cat) + Mono + Gran3cat + 
                     Bcell + CD8T + CD4T, 
                   data = data)


# check assumptions
test_ph <- cox.zph(cox_model)
print(test_ph)

print(ifelse(test_ph$table[,3] < 0.05, "yes", "no")[test_ph$table[,3] < 0.05])
# if p-value is significant, stratify by that variable
# p-values < 0.05 for Gran, NK, and BMI

# Summarize the model
cox_summary <- summary(cox_model)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "OR", "SE",
                            "z", "P_val")

write.csv(coefficients, file = "v2_cox_results_fdr.csv")

ggforest(cox_model, data = data)

####### Generate Kaplan-Meier curve #########

# Create survival object
surv_obj <- Surv(data$days_to_event, data$CHD == 1)

# Fit Kaplan-Meier survival curves
km_fit <- survfit(surv_obj ~ stress2cat, data = data)

event_plot <- ggsurvplot(km_fit, 
                         data = data, 
                         pval = TRUE, 
                         conf.int = TRUE,
                         fun = "event",  # This inverts to show event probability
                         legend.labs = c("Low", "High"),
                         legend.title = "Stress Level",
                         title = "CHD Cumulative Incidence by Stress level",
                         xlab = "Time to Event in Days",
                         ylab = "CHD Probability")$plot

event_plot

#######################################
###### Threshold 2 (Bonferroni) #######
#######################################

data_unfiltered <- read.csv("chd_data_v2.csv", header = TRUE)

data <- data_unfiltered %>% drop_na(MRS2cat_bonf) %>% drop_na(CHD) %>% drop_na(days_to_event) %>%
  drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore)
dim(data)
#[1] 1202   43

# Ensure MRS_med_split is a factor with proper labels
data$MRS2cat_bonf <- factor(data$MRS2cat_bonf, 
                       levels = c(1, 2), 
                       labels = c("Low MRS", "High MRS"))

data$BMI3cat <- ntile(data$BMI, 3)
data$NK3cat <- ntile(data$NK, 3)
data$Gran3cat <- ntile(data$Gran, 3)

######## Fit Cox proportional hazards regression model #########
cox_model <- coxph(Surv(days_to_event, CHD) ~ MRS2cat_bonf + sex + age + 
                     currentSmoker + strata(BMI3cat) + Income + edu3cat + alc + 
                     marital2cat + nbSESanascore + 
                     PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                     PC8 + PC9 + PC10 + strata(NK3cat) + Mono + Gran3cat + 
                     Bcell + CD8T + CD4T, 
                   data = data)

# check assumptions
test_ph <- cox.zph(cox_model)
print(test_ph)

print(ifelse(test_ph$table[,3] < 0.05, "yes", "no")[test_ph$table[,3] < 0.05])
# BMI, NK, and Gran are significant - stratify into tertiles

# Summarize the model
cox_summary <- summary(cox_model)

# Summarize the model
coefficients <- summary(cox_model)$coefficients
colnames(coefficients) <- c("Estimate", "OR", "SE",
                            "z", "P_val")

write.csv(coefficients, file = "v2_cox_results_bonferroni.csv")

# In WHI, MRS is negatively associated. Higher in stress group, 
# closer to 0, less negative, visualize association with stress 

# check assumptions
test_ph <- cox.zph(cox_model)
print(test_ph)
# if p-value is significant, you need to stratify by that variable
# BMI was significant when treated as a continuous variable,
# so I split BMI into 3 categories
# BMI did not have a significant p-value after splitting into tertiles

ggforest(cox_model, data = data)

####### Generate Kaplan-Meier curve #########

# Create survival object
surv_obj <- Surv(data$days_to_event, data$CHD == 1)

# Fit Kaplan-Meier survival curves
km_fit <- survfit(surv_obj ~ stress2cat, data = data)

event_plot <- ggsurvplot(km_fit, 
                         data = data, 
                         pval = TRUE, 
                         conf.int = TRUE,
                         fun = "event",  # This inverts to show event probability
                         legend.labs = c("Low", "High"),
                         legend.title = "Stress Level",
                         title = "CHD Cumulative Incidence by Stress level",
                         xlab = "Time to Event in Days",
                         ylab = "CHD Probability")$plot

event_plot
