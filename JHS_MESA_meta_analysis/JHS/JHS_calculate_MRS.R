##################################################
########### JHS EPICv1 MRS calculation ###########
##################################################

library(dplyr)
library(tidyverse)

path = "/work/users/h/e/hemilla"
path2 = "/work/users/h/e/hemilla/WHI_JHS_MRS"
setwd(path)

# get beta values for JHS data
load("beta.noXY.RData")

# get coefficients to calculate methylation risk scores
setwd(path2)
coef <- read.csv("meta_coef_sig_annotated_SLEs.csv", header = TRUE,
                 row.names = 1)

####### THRESHOLD 1 #######
# make estimates & beta values compatible for linear combination
beta <- beta.noXY[rownames(beta.noXY) %in% coef$IlmnID,]
est <- coef[coef$IlmnID %in% rownames(beta),]
dim(est) #770 x 14

est <- est[order(match(est$IlmnID,rownames(beta))),]
identical(est$IlmnID, rownames(beta))

estimates <- as.matrix(est$beta)
beta_matrix <- t(as.matrix(beta))

# Calculate the score using matrix multiplication
MRS <-beta_matrix %*% estimates

# Add sample names to the score
names(MRS) <- rownames(beta_matrix)

summary(MRS)

###### CALCULATE MRS FOR BONFERRONI-CORRECTED SIGNIFICANT CPG SITES #########

coef_bonf <- filter(coef, bonferroni < 0.05)
dim(coef_bonf)
#[1] 13 14

# make estimates & beta values compatible for linear combination
beta <- beta.noXY[rownames(beta.noXY) %in% coef_bonf$IlmnID,]
est <- coef_bonf[coef_bonf$IlmnID %in% rownames(beta),]

est <- est[order(match(est$IlmnID,rownames(beta))),]
identical(est$IlmnID, rownames(beta))

estimates <- as.matrix(est$beta)
beta_matrix <- t(as.matrix(beta))

# Calculate the score using matrix multiplication
MRS_bonferroni <- beta_matrix %*% estimates

# Add sample names to the score
names(MRS) <- rownames(beta_matrix)

# make data frame with MRS values for all thresholds

df <- data.frame(MRS = MRS, MRS_FDR = MRS_FDR, 
                    MRS_bonferroni = MRS_bonferroni, 
                    row.names = rownames(beta_matrix))

# write file with MRS values
write.csv(df, file = "MRS_JHS_EPICV1.csv")

##################################################
########### JHS EPICv2 MRS calculation ###########
##################################################

######## CALCULATE MRS FOR 840 CPG SITES FROM WHI ###########
# coefficients related to stress calculated for WHI

library(dplyr)
library(tidyverse)

path = "/work/users/h/e/hemilla/WHI_JHS_MRS/EPICv2"
setwd(path)

# remove samples that overlap v1 data
# see lines 1-21 of "calculate_MRS_040724.R" to determine which samples
# were used in the v1 analysis

# get correspondence between Sample_Name and SUBJECT_ID
#pheno <- read.csv("pheno_020524.csv", header = TRUE)

#samples_v1 <- data.frame(Sample_Name = colnames(beta) %>% as.integer())
# 1709 participants
#sample <- merge(samples_v1, pheno, by = "Sample_Name")
#write.csv(sample$SUBJECT_ID, "samples_used_subjid_MRS.csv", row.names = FALSE)

sample <- read.csv("samples_used_subjid_MRS.csv", header = TRUE)
dim(sample)
#1691

#Use this beta.v2.subjid.RData prepared in the file "FVIII_v2_data_prep.R"
#beta.v2.subjid.RData uses subject_ID instead of TOEID

# get beta values for JHS data (v2)
load("beta.v2.subjid.RData")

overlapping_samples <- colnames(beta)[colnames(beta) %in% sample$x]
length(overlapping_samples)
#[1] 67

beta1 <- beta[, !(colnames(beta) %in% overlapping_samples)] # remove overlapping samples
dim(beta)
#[1] 860960   1687
dim(beta1)
#[1] 860960   1620

# get coefficients to calculate methylation risk scores
coef <- read.csv("meta_coef_sig_annotated_SLEs.csv", header = TRUE,
                 row.names = 1)

####### THRESHOLD 1 #######
# make estimates & beta values compatible for linear combination
beta1 <- beta1[rownames(beta1) %in% coef$IlmnID,]
est <- coef[coef$IlmnID %in% rownames(beta1),]
dim(est) #773  14

est <- est[order(match(est$IlmnID,rownames(beta1))),]
identical(est$IlmnID, rownames(beta1))

estimates <- as.matrix(est$beta)
beta_matrix <- t(as.matrix(beta1))

# Calculate the score using matrix multiplication
MRS <-beta_matrix %*% estimates

# Add sample names to the score
names(MRS) <- rownames(beta_matrix)

summary(MRS)

###### CALCULATE MRS FOR BONFERRONI-CORRECTED SIGNIFICANT CPG SITES #########

coef_bonf <- filter(coef, bonferroni < 0.05)
dim(coef_bonf)
#[1] 13 14

# make estimates & beta values compatible for linear combination
beta_bonf <- beta1[rownames(beta1) %in% coef_bonf$IlmnID,]
dim(beta_bonf) #[1]   12 1620
est <- coef_bonf[coef_bonf$IlmnID %in% rownames(beta_bonf),]
dim(est) #[1] 12 14

est <- est[order(match(est$IlmnID,rownames(beta_bonf))),]
identical(est$IlmnID, rownames(beta_bonf))

estimates <- as.matrix(est$beta)
beta_matrix <- t(as.matrix(beta_bonf))

# Calculate the score using matrix multiplication
MRS_bonferroni <- beta_matrix %*% estimates

# Add sample names to the score
names(MRS_bonferroni) <- rownames(beta_matrix)
summary(MRS_bonferroni)

# make data frame with MRS values for all thresholds

df <- data.frame(SUBJECT_ID = rownames(beta_matrix), 
                 MRS = MRS, MRS_FDR = MRS_FDR, 
                 MRS_bonferroni = MRS_bonferroni)

# write file with MRS values
write.csv(df, file = "MRS_JHS_EPICV2.csv", row.names = FALSE)
