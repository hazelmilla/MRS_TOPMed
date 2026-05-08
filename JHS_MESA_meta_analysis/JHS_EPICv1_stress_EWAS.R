# Run analysis for individual CpGs (beyond MRS) & stress only. 
# How many CpGs (of 841) individually replicate in association with stress in the JHS cohort
library(dplyr)
library(tidyverse)

path = "/work/users/h/e/hemilla"
setwd(path)

# prepare data
data <- read.csv("prepared_data.csv", 
                     header = TRUE)
load("beta.noXY.RData")

path2 = "/work/users/h/e/hemilla/WHI_JHS_MRS"
setwd(path2)
coef <- read.csv("meta_coef_sig_annotated_SLEs.csv", header = TRUE,
                 row.names = 1)

###### All participants ########
# see if individual CpGs are associated with stress
# sex, age, smoking, alc, BMI, marriage, SES, PC1-10, Houseman cell types
data <- data %>% drop_na(stress2cat) %>% drop_na(sex) %>% drop_na(age) %>%
  drop_na(currentSmoker) %>% drop_na(BMI) %>% drop_na(Income) %>% 
  drop_na(edu3cat) %>% drop_na(alc) %>% drop_na(marital2cat) %>% 
  drop_na(nbSESanascore) %>%
  drop_na(PC1) %>% drop_na(PC2) %>% drop_na(PC3) %>% 
  drop_na(PC4) %>% drop_na(PC5) %>% drop_na(PC6) %>% drop_na(PC7) %>% 
  drop_na(PC8) %>% drop_na(PC9) %>% drop_na(PC10) %>% drop_na(NK) %>% 
  drop_na(Mono) %>% drop_na(Gran) %>% drop_na(Bcell) %>% drop_na(CD8T) %>%
  drop_na(CD4T)
dim(data)
#[1] 1417   37

beta <- beta.noXY[rownames(beta.noXY) %in% coef$IlmnID, colnames(beta.noXY) %in% data$Sample_Name]
beta <- beta[, order(match(colnames(beta), data$Sample_Name))]
identical((colnames(beta1) %>% as.integer()), data1$Sample_Name)
dim(beta) #[1]  770 1417

# Create a design matrix
design <- model.matrix(~ stress2cat + sex + age + currentSmoker + 
                         BMI + Income + edu3cat + alc + marital2cat +
                         nbSESanascore +
                         PC1 + PC2 +  PC3 + PC4 + PC5 + PC6 + PC7 + 
                         PC8 + PC9 + PC10 + NK + Mono + Gran + 
                         Bcell + CD8T + CD4T, data=data) 
design <- design %>% as.data.frame()

results <- matrix(data = NA, ncol = 4, nrow = length(beta[,1]))

for (i in 1:length(beta[,1])){
  
fit <- lm(beta[i,] ~ stress2cat + sexMale + age + currentSmoker + 
            BMI + Income + edu3cat + alc + marital2catNot_married +
            nbSESanascore + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
            PC7 + PC8 + PC9 + PC10 + NK + Mono + Gran + Bcell + 
            CD8T + CD4T, data=design)
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

rownames(results) <- rownames(beta)
colnames(results) <- c("Estimate", "Std.Error", "t_value",
                       "P_value")

results_df <- results %>% as.data.frame()

results_df1 <- results_df %>%
  mutate(bonferroni = p.adjust(P_value, method="bonferroni"),
         hochberg = p.adjust(P_value, method="hochberg"))
sum(results_df1$bonferroni < 0.05) #0
sum(results_df1$hochberg < 0.05) #0
sum(results_df1$P_value < 0.05) #82

write.csv(results_df1, file = "stress_WHI_CpG_lm_results_Apr25.csv")

# compare significant sites coefficient directions in WHI & JHS

library(tidyverse)
library(dplyr)

path2 = "/work/users/h/e/hemilla/WHI_JHS_MRS"
setwd(path2)
coef <- read.csv("meta_coef_sig_annotated_SLEs.csv", header = TRUE,
                 row.names = 1)

path = "/work/users/h/e/hemilla/WHI_JHS_MRS/stress_CpG_sites_WHI-JHS"
setwd(path)

jhs <- read.csv("stress_WHI_CpG_lm_results_Apr25.csv", header = TRUE)
whi <- coef

(whi$beta)
sign_jhs <- filter(jhs, P_value < 0.05)

# Filter to keep only matching IlmnIDs
whi1 <- whi %>% filter(IlmnID %in% sign_jhs$X)
# Reorder rows to match sign_jhs$X
whi1 <- whi1 %>% arrange(match(IlmnID, sign_jhs$X))
# Verify order
identical(sign_jhs$X, whi1$IlmnID)
# [1] TRUE

sign_matches <- sign(sign_jhs$Estimate) == sign(whi1$beta)

sum(sign_matches)
#[1] 77 (of 82)

# Determine direction of significance in female participants
jhs <- read.csv("stress_WHI_CpG_lm_female_Apr25.csv", header = TRUE)
whi <- coef

(whi$beta)
sign_jhs <- filter(jhs, P_value < 0.05)

# Filter to keep only matching IlmnIDs
whi1 <- whi %>% filter(IlmnID %in% sign_jhs$X)
# Reorder rows to match sign_jhs$X
whi1 <- whi1 %>% arrange(match(IlmnID, sign_jhs$X))
# Verify order
identical(sign_jhs$X, whi1$IlmnID)
# [1] TRUE

sign_matches <- sign(sign_jhs$Estimate) == sign(whi1$beta)

sum(sign_matches) 
# [1] 154 (out of 155)

# check for sites that are significant after Bonferroni correction
bonf_jhs <- filter(jhs, bonferroni < 0.05)
# Filter to keep only matching IlmnIDs
whi1 <- whi %>% filter(IlmnID %in% bonf_jhs$X)
# Reorder rows to match bonf_jhs$X
whi1 <- whi1 %>% arrange(match(IlmnID, bonf_jhs$X))
# Verify order
identical(bonf_jhs$X, whi1$IlmnID)
# [1] TRUE

sign_matches <- sign(bonf_jhs$Estimate) == sign(whi1$beta)
sum(sign_matches) 
# [1] 2 (of 2)

# Determine direction of significance in male participants
jhs <- read.csv("stress_WHI_CpG_lm_male_Apr25.csv", header = TRUE)
whi <- coef

(whi$beta)
sign_jhs <- filter(jhs, P_value < 0.05)

# Filter to keep only matching IlmnIDs
whi1 <- whi %>% filter(IlmnID %in% sign_jhs$X)
# Reorder rows to match sign_jhs$X
whi1 <- whi1 %>% arrange(match(IlmnID, sign_jhs$X))
# Verify order
identical(sign_jhs$X, whi1$IlmnID)
# [1] TRUE

sign_matches <- sign(sign_jhs$Estimate) == sign(whi1$beta)

sum(sign_matches) 
# [1] 9 (of 19)
