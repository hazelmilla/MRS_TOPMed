### Calculate MRS for 3 thresholds ###
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

# make data frame with MRS values for all 3 thresholds

df <- data.frame(MRS = MRS, MRS_FDR = MRS_FDR, 
                    MRS_bonferroni = MRS_bonferroni, 
                    row.names = rownames(beta_matrix))

# write file with MRS values
write.csv(df, file = "MRS_JHS_EPICV1.csv")
