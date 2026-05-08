# Run linear regression between individual CpG sites and stress score

library(dplyr)
library(tidyverse)

path = "/work/users/h/e/hemilla/WHI_JHS_MRS/EPICv2"
setwd(path)

# prepare data
data <- read.csv("prepared_data_v2.csv", 
                 header = TRUE)
coef <- read.csv("meta_coef_sig_annotated_SLEs.csv", header = TRUE,
                 row.names = 1)
load("beta.v2.subjid.RData")

# remove overlapping samples
sample <- read.csv("samples_used_subjid_MRS.csv", header = TRUE)
dim(sample)
#[1] 1691    1
overlapping_samples <- colnames(beta)[colnames(beta) %in% sample$x]
length(overlapping_samples)
#[1] 67
# remove overlapping samples
beta1 <- beta[, !(colnames(beta) %in% overlapping_samples)]
dim(beta)
#[1] 860960   1687
dim(beta1)
#[1] 860960   1620

path2 = "/work/users/h/e/hemilla/WHI_JHS_MRS/EPICv2/indiv_CpG_stress"
setwd(path2)

# see if individual CpGs are associated with stress
# covariates: sex, age, smoking, alc, BMI, marriage, SES, PC1-10, Houseman cell types
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
#[1] 1312   36

beta <- beta1[rownames(beta1) %in% coef$IlmnID , colnames(beta1) %in% data$SUBJECT_ID]
beta <- beta[, order(match(colnames(beta), data$SUBJECT_ID))]
data <- data[data$SUBJECT_ID %in% colnames(beta),]
identical((colnames(beta) %>% as.integer()), data$SUBJECT_ID)
dim(beta) #[1]  773 1312

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
sum(results_df1$P_value < 0.05) #17

write.csv(results_df1, file = "stress_WHI_CpG_lm_results_v2_Apr25.csv")
