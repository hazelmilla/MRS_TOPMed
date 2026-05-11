#####################################################
############# THRESHOLD 1 - FDR ##############
#####################################################

# Run analysis for individual CpGs (beyond MRS) & stress only. 
# How many CpGs (of 841) individually replicate in terms of 
# their association with stress in the JHS cohort?
library(dplyr)
library(tidyverse)

path = "/work/users/h/e/hemilla/WHI_MESA_MRS"
setwd(path)

data_unfiltered <- read.csv("Input/prepared_data_MESA.csv", header = TRUE)

# get beta values
load("Input/MESA_methyl_visit_1_cleaned.RData")
dim(methylation.vis1)

# get only CpG beta values
beta <- methylation.vis1[, grep("^cg", colnames(methylation.vis1))]
beta <- t(beta) # transpose beta so rows are CpGs and cols are Sample ID
# get all other info 
dat1 <- methylation.vis1[, -(grep("^cg", colnames(methylation.vis1)))]
dat1 <- dat1 %>% as.data.frame()

coef <- read.csv("Input/meta_coef_sig_annotated_SLEs.csv", header = TRUE,
                 row.names = 1)

var <- c("stress2cat")
cov <- c("gender1", "AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
         "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
         "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", "Sample_Plate", 
         "Sample_Well", "Sentrix_ID", "NK", "Gran", "Mono",
         "Bcell", "CD8T", "CD4T")
length(cov)
#29

data <- data_unfiltered[complete.cases(data_unfiltered[ , cov]) & complete.cases(data_unfiltered[ , var]), ]
dim(data)
#[1] 698  93

beta1 <- beta[rownames(beta) %in% coef$IlmnID, colnames(beta) %in% data$TOEID]
beta1 <- beta1[, order(match(colnames(beta1), data$TOEID))]
identical(colnames(beta1), data$TOEID) # TRUE
dim(beta) # [1] 768840    870
dim(beta1) # [1] 784 600

# Run linear regression

results <- matrix(data = NA, ncol = 4, nrow = length(beta1[,1]))

for (i in 1:length(beta1[,1])){
  
  data$beta_val <- beta1[i, ]
  
  formula_str <- paste("beta_val ~", "stress2cat", "+", paste(cov, collapse = " + "))
  model_formula <- as.formula(formula_str)
  
  fit <- lm(model_formula, data=data)
  sumfit <- summary(fit)
  
  coefficients <- sumfit$coefficients
  colnames(coefficients) <- c("Estimate", "Std.Error", 
                              "t_value", "P_value")
  
  vect <- c(coefficients["stress2cat","Estimate"], 
            coefficients["stress2cat","Std.Error"], 
            coefficients["stress2cat","t_value"],
            coefficients["stress2cat","P_value"])
  results[i,] <- vect
}


rownames(results) <- rownames(beta1)
colnames(results) <- c("Estimate", "Std.Error", "t_value",
                       "P_value")

results_df <- results %>% as.data.frame()

results_df1 <- results_df %>%
  mutate(bonferroni = p.adjust(P_value, method="bonferroni"),
         hochberg = p.adjust(P_value, method="hochberg"))
sum(results_df1$bonferroni < 0.05) #0
sum(results_df1$hochberg < 0.05) #0
sum(results_df1$P_value < 0.05) #24

mesa <- results_df1
whi <- coef
sign_mesa <- filter(mesa, P_value < 0.05)

# Filter to keep only matching IlmnIDs
whi1 <- whi %>% filter(IlmnID %in% rownames(sign_mesa))
# Reorder rows to match sign_jhs$X
whi1 <- whi1 %>% arrange(match(IlmnID, rownames(sign_mesa)))
# Verify order
identical(rownames(sign_mesa), whi1$IlmnID)
# [1] TRUE

sign_matches <- sign(sign_mesa$Estimate) == sign(whi1$beta)

sum(sign_matches) # 10 matches

write.csv(results_df1, file = "Output/MESA_stress_WHI_CpGsites_lm.csv")

#####################################################
############# THRESHOLD 2 - BONFERRONI ##############
#####################################################

library(dplyr)
library(tidyverse)

path = "/work/users/h/e/hemilla/WHI_MESA_MRS"
setwd(path)

data_unfiltered <- read.csv("Input/prepared_data_MESA.csv", header = TRUE)

# get beta values
load("Input/MESA_methyl_visit_1_cleaned.RData")
dim(methylation.vis1) #[1]    870 771361

# get only CpG beta values
beta <- methylation.vis1[, grep("^cg", colnames(methylation.vis1))]
beta <- t(beta) # transpose beta so rows are CpGs and cols are Sample ID

coef <- read.csv("Input/meta_coef_sig_annotated_SLEs.csv", header = TRUE,
                 row.names = 1)

coef_bonf <- filter(coef, coef$bonferroni < 0.05)
dim(coef_bonf) # [1] 13 14

var <- c("stress2cat")
cov <- c("gender1", "AGE_AT_COLLECTION", "bmi1c", "marital2cat", 
         "income4cat", "edu4cat", "race1c", "F1_PC2_1", "curalc1", "cursmk1", "pc1", "pc2", 
         "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", "Sample_Plate", 
         "Sample_Well", "Sentrix_ID", "NK", "Gran", "Mono",
         "Bcell", "CD8T", "CD4T")

data <- data_unfiltered[complete.cases(data_unfiltered[ , cov]) & complete.cases(data_unfiltered[ , var]), ]
dim(data)
#[1] 698  94

beta1 <- beta[rownames(beta) %in% coef_bonf$IlmnID, colnames(beta) %in% data$TOEID]
beta1 <- beta1[, order(match(colnames(beta1), data$TOEID))]
identical(colnames(beta1), data$TOEID) # TRUE
dim(beta) # [1] 768840    870
dim(beta1) # [1] 12 698

# Run linear regression
results <- matrix(data = NA, ncol = 4, nrow = length(beta1[,1]))

for (i in 1:length(beta1[,1])){
  
  data$beta_val <- beta1[i, ]
  
  formula_str <- paste("beta_val ~", "stress2cat", "+", paste(cov, collapse = " + "))
  model_formula <- as.formula(formula_str)
  
  fit <- lm(model_formula, data=data)
  sumfit <- summary(fit)
  
  coefficients <- sumfit$coefficients
  colnames(coefficients) <- c("Estimate", "Std.Error", 
                              "t_value", "P_value")
  
  vect <- c(coefficients["stress2cat","Estimate"], 
            coefficients["stress2cat","Std.Error"], 
            coefficients["stress2cat","t_value"],
            coefficients["stress2cat","P_value"])
  results[i,] <- vect
}


rownames(results) <- rownames(beta1)
colnames(results) <- c("Estimate", "Std.Error", "t_value",
                       "P_value")

results_df <- results %>% as.data.frame()

results_df1 <- results_df %>%
  mutate(bonferroni = p.adjust(P_value, method="bonferroni"),
         hochberg = p.adjust(P_value, method="hochberg"))
sum(results_df1$bonferroni < 0.05) #0
sum(results_df1$hochberg < 0.05) #0
sum(results_df1$P_value < 0.05) #0

write.csv(results_df1, file = "Output/MESA_stress_WHI_CpGsites_lm_Bonferroni.csv")
