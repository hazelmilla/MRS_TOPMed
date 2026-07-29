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

# See lines 1-37 to get list of common TNF sites across cohorts
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
