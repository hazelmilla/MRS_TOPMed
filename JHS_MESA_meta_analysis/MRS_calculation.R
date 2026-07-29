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
dim(beta_matrix) # 1709   12
# 12 of the 13 Bonferroni sites are present

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

######## Calculate MRS for TNF sites only ########
library(dplyr)
library(tidyverse)

path = "/work/users/h/e/hemilla/JHS_stress_CpG/Input_files"
setwd(path)

load("beta.noXY.RData")
load("beta.v2.subjid.RData")
beta_v2 <- beta

coef <- read.csv("meta_coef_sig_annotated_SLEs.csv", header = TRUE,
                 row.names = 1)

setwd("/work/users/h/e/hemilla/WHI_MESA_MRS/Input")
load("MESA_methyl_visit_1_cleaned.RData")
beta_mesa <- methylation.vis1[, grep("^cg", colnames(methylation.vis1))]
beta_mesa <- t(beta_mesa)

common_probes <- Reduce(intersect, list(rownames(beta.noXY), rownames(beta_v2), rownames(beta_mesa)))
length(common_probes) #[1] 659563

tnf_all <- c("cg13868520", "cg10594075", "cg13203480", 
             "cg18064152", "cg02137984", "cg04472685", 
             "cg05952498", "cg09637172", "cg17755321", 
             "cg20477259", "cg24452282", "cg26786341", 
             "cg27531490")
tnf_no_receptors <- c("cg02137984", "cg04472685", 
                      "cg05952498", "cg09637172", "cg17755321", 
                      "cg20477259", "cg24452282", "cg26786341", 
                      "cg27531490")

tnf_noR_common <- common_probes[common_probes %in% tnf_no_receptors]
length(tnf_noR_common)

beta <- beta.noXY[rownames(beta.noXY) %in% tnf_noR_common,]
est <- coef[coef$IlmnID %in% rownames(beta),]
est <- est[order(match(est$IlmnID,rownames(beta))),]

dim(beta) #[1]    4 1709
dim(est) #[1]  4 14

identical(est$IlmnID, rownames(beta)) # TRUE

estimates <- as.matrix(est$beta)
beta_matrix <- t(as.matrix(beta))

# Calculate the score using matrix multiplication
MRS <-beta_matrix %*% estimates
names(MRS) <- rownames(beta_matrix)
summary(MRS)

# Calculate TNF average
TNF_avg <- colMeans(beta)

df <- data.frame(TNF_MRS = MRS, TNF_avg = TNF_avg)

setwd("/work/users/h/e/hemilla/JHS_stress_CpG/Output_files")
write.csv(df, "TNF_sites_MRS_JHSv1.csv", row.names = TRUE)

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

######## Calculate MRS for TNF sites only ########
# remove samples that overlap v1 data
# see lines 1-21 of "calculate_MRS_040724.R" to determine which samples
# were used in the v1 analysis

#sample <- read.csv("Input_files/samples_used_subjid_MRS.csv", header = TRUE)
#dim(sample)
#1691

#Use this beta.v2.subjid.RData prepared in the file "FVIII_v2_data_prep.R"
#beta.v2.subjid.RData uses subject_ID instead of TOEID

# get beta values for JHS data (v2)
# load("Input_files/beta.v2.subjid.RData")

#overlapping_samples <- colnames(beta)[colnames(beta) %in% sample$x]
#length(overlapping_samples)
#[1] 67
# remove overlapping samples
#beta1 <- beta[, !(colnames(beta) %in% overlapping_samples)]
#dim(beta)
#[1] 860960   1687
#dim(beta1)
#[1] 860960   1620

# make estimates & beta values compatible for linear combination
beta1 <- beta_v2[rownames(beta_v2) %in% tnf_noR_common,]
est <- coef[coef$IlmnID %in% rownames(beta1),]
dim(est) #4 14
dim(beta1) #4 1687

est <- est[order(match(est$IlmnID,rownames(beta1))),]
identical(est$IlmnID, rownames(beta1))

estimates <- as.matrix(est$beta)
beta_matrix <- t(as.matrix(beta1))

# Calculate the score using matrix multiplication
MRS <-beta_matrix %*% estimates

# Add sample names to the score
names(MRS) <- rownames(beta_matrix)
summary(MRS)

# Calculate TNF average
TNF_avg <- colMeans(beta1)

df <- data.frame(TNF_MRS = MRS, TNF_avg = TNF_avg)
write.csv(df, "TNF_sites_MRS_JHSv2.csv", row.names = TRUE)

###############################################
########### Calculate MRS for MESA ############
###############################################
# MESA calculate MRS score #

# filter for common TNF sites across MESA and JHS

library(tidyverse)
library(dplyr)
library(data.table)

path = "/work/users/h/e/hemilla/WHI_MESA_MRS"
setwd(path)

# get beta values
load("Input/MESA_methyl_visit_1_cleaned.RData")
dim(methylation.vis1)

# get only CpG beta values
beta <- methylation.vis1[, grep("^cg", colnames(methylation.vis1))]
# get all other info 
dat <- methylation.vis1[, -(grep("^cg", colnames(methylation.vis1)))]
dat <- dat %>% as.data.frame()

# get coefficients to calculate methylation risk scores
coef <- read.csv("Input/meta_coef_sig_annotated_SLEs.csv", header = TRUE,
                 row.names = 1)

####### THRESHOLD 1 (FDR) #######
# check if you need to transpose beta matrix
dim(beta)
beta <- t(beta)
rownames(beta)
colnames(beta)

# make estimates & beta values compatible for linear combination
beta1 <- beta[rownames(beta) %in% coef$IlmnID,]
est <- coef[coef$IlmnID %in% rownames(beta1),]
dim(est) #[1] 784  14

est <- est[order(match(est$IlmnID,rownames(beta1))),]
identical(est$IlmnID, rownames(beta1))

estimates <- as.matrix(est$beta)
beta_matrix <- t(as.matrix(beta1))

# Calculate the score using matrix multiplication
MRS <- beta_matrix %*% estimates

# Add sample names to the score
names(MRS) <- rownames(beta_matrix)

summary(MRS)

###### THRESHOLD 2 (BONFERRONI) ######
coef_bonf <- filter(coef, bonferroni < 0.05)
dim(coef_bonf)
#[1] 13 14

# make estimates & beta values compatible for linear combination
beta2 <- beta[rownames(beta) %in% coef_bonf$IlmnID,]
est <- coef_bonf[coef_bonf$IlmnID %in% rownames(beta2),]

est <- est[order(match(est$IlmnID,rownames(beta2))),]
identical(est$IlmnID, rownames(beta2))
dim(beta2) #[1] 12 14
dim(est) #[1]  12 870

estimates <- as.matrix(est$beta)
beta_matrix <- t(as.matrix(beta2))

# Calculate the score using matrix multiplication
MRS_bonf <- beta_matrix %*% estimates

# Add sample names to the score
names(MRS_bonf) <- rownames(beta_matrix)
summary(MRS_bonf)

# Create dataframe with thresholds
df <- data.frame(MRS = MRS, MRS_bonf = MRS_bonf, 
                 row.names = rownames(beta_matrix))

# write file with MRS values
write.csv(df, file = "Input/MRS_MESA.csv")


####### Calculate MRS for TNF receptor sites #########

#library(dplyr)
#library(tidyverse)

#path = "/work/users/h/e/hemilla/WHI_MESA_MRS"
#setwd(path)

# get beta values
#load("Input/MESA_methyl_visit_1_cleaned.RData")
#dim(methylation.vis1)
# get only CpG beta values
#beta <- methylation.vis1[, grep("^cg", colnames(methylation.vis1))]
#beta <- t(beta)

#coef <- read.csv("Input/meta_coef_sig_annotated_SLEs.csv", header = TRUE,
#                 row.names = 1)

# See lines 95-106 to get list of common TNF sites across cohorts
tnf_noR_common

beta <- beta_mesa[rownames(beta_mesa) %in% tnf_noR_common,]
est <- coef[coef$IlmnID %in% rownames(beta),]
est <- est[order(match(est$IlmnID,rownames(beta))),]

dim(beta) #[1]   4 870
dim(est) #[1]  4 14

identical(est$IlmnID, rownames(beta))

estimates <- as.matrix(est$beta)
beta_matrix <- t(as.matrix(beta))

# Calculate the score using matrix multiplication
MRS <-beta_matrix %*% estimates
names(MRS) <- rownames(beta_matrix)
summary(MRS)

# Calculate TNF average
TNF_avg <- colMeans(beta)

df <- data.frame(MRS_TNF = MRS, TNF_avg = TNF_avg)

setwd("/work/users/h/e/hemilla/WHI_MESA_MRS/Output")
write.csv(df, "TNF_sites_MRS_MESA.csv", 
          row.names = TRUE)
