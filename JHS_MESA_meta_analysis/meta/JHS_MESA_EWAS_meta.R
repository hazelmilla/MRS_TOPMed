##### JHS MESA meta-analysis ######

library(metafor)
library(tidyverse)

setwd("/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/Meta_analysis")

###### Individual CpG sites and stress #######
v1 <- read.csv("v1_stress_WHI_CpG_lm_results_Apr25.csv", header = TRUE)
v2 <- read.csv("v2_stress_WHI_CpG_lm_results_Apr25.csv", header = TRUE)


setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis")
mesa <- read.csv("Input/MESA_stress_WHI_CpGsites_lm.csv", header = TRUE)

setwd("/Users/hazelmilla/Desktop/Zannas_lab/JHS_MESA_meta_analysis")

v1$study <- ("JHSv1")
v2$study <- ("JHSv2")
mesa$study <- ("MESA")

colnames(mesa) <- c("CpG", "BETA", "SE", "t_value", "P_VAL", "bonferroni",
                    "hochberg", "study")

# Step 1: Find common row names
common_rows <- Reduce(intersect, list(mesa$CpG, v1$CpG, v2$CpG))
length(common_rows)
#[1] 731

# Step 2: Subset each data frame by those row names
mesa_common <- mesa[mesa$CpG %in% common_rows, , drop = FALSE]
v1_common <- v1[v1$CpG %in% common_rows, , drop = FALSE]
v2_common <- v2[v2$CpG %in% common_rows, , drop = FALSE]

mesa_dat <- mesa_common[match(common_rows, mesa_common$CpG),]
v1_dat <- v1_common[match(common_rows, v1_common$CpG),]
v2_dat <- v2_common[match(common_rows, v2_common$CpG),]

identical(v1_dat$CpG, v2_dat$CpG) 
identical(v1_dat$CpG, mesa_dat$CpG)

sitebind <- NULL
run_meta <- NULL
run_coef <- NULL
meta_coef <- NULL
i <- 1

for (i in (1:length(v1_dat$CpG))){
  sitebind <- rbind(v1_dat[i,], v2_dat[i,], mesa_dat[i,])
  run_meta <- rma(yi = sitebind$BETA, sei = sitebind$SE, slab = sitebind$study, method = "FE")
  run_coef <- c(run_meta$b, run_meta$se, run_meta$zval, run_meta$QEp, run_meta$pval, run_meta$ci.lb, run_meta$ci.ub)
  meta_coef <- rbind(meta_coef, run_coef)
  print(i)
  flush.console()
}

dim(meta_coef)
class(meta_coef)
meta_coef_df <- as.data.frame(meta_coef)
colnames(meta_coef_df) <- c("beta", "SE", "Z", "p_het", "p_meta", "lower CI", "upper CI")
rownames(meta_coef_df) <- v1_dat$CpG

meta_coef_df$q_meta <- p.adjust(meta_coef_df$p_meta, method = "fdr")
meta_coef_df$bonferroni <- p.adjust(meta_coef_df$p_meta, method = "bonferroni")

sum(meta_coef_df$p_het < 0.05) #29
sum(meta_coef_df$p_meta < 0.05) # 57
sum(meta_coef_df$q_meta < 0.05) # 0
sum(meta_coef_df$bonferroni < 0.05) # 0

write.csv(meta_coef_df, file = "Output/CpG_stress_meta_MESA_JHS_all_common_probes_FDR.csv")


# run CpGs with significant heterogeneity for random effects
sig_hg <- meta_coef_df %>% filter(meta_coef_df$p_het < 0.05)
sig_hg <- rownames(sig_hg)

mesa_cp <- mesa[mesa$CpG %in% sig_hg,]
v1_cp <- v1[v1$CpG %in% sig_hg,]
v2_cp <- v2[v2$CpG %in% sig_hg,]

mesa_dat <- mesa_cp[match(sig_hg, mesa_cp$CpG),]
v1_dat <- v1_cp[match(sig_hg, v1_cp$CpG),]
v2_dat <- v2_cp[match(sig_hg, v2_cp$CpG),]
identical(v1_dat$CpG, v2_dat$CpG)
identical(v1_dat$CpG, mesa_dat$CpG)

sitebind <- NULL
run_meta <- NULL
run_coef <- NULL
meta_coef <- NULL
i <- 1

for (i in (1:length(v1_dat$CpG))){
  sitebind <- rbind(v1_dat[i,], v2_dat[i,], mesa_dat[i,])
  run_meta <- rma(yi = sitebind$BETA, sei = sitebind$SE, slab = sitebind$study, method = "REML")
  run_coef <- c(run_meta$b, run_meta$se, run_meta$zval, run_meta$QEp, run_meta$pval, run_meta$ci.lb, run_meta$ci.ub)
  meta_coef <- rbind(meta_coef, run_coef)
  print(i)
  flush.console()
}

dim(meta_coef)
class(meta_coef)
meta_coef_df <- as.data.frame(meta_coef)
colnames(meta_coef_df) <- c("beta", "SE", "Z", "p_het", "p_meta", "lower CI", "upper CI")
rownames(meta_coef_df) <- v1_dat$CpG

meta_coef_df$q_meta <- p.adjust(meta_coef_df$p_meta, method = "fdr")
meta_coef_df$bonferroni <- p.adjust(meta_coef_df$p_meta, method = "bonferroni")

sum(meta_coef_df$p_het < 0.05) #29
sum(meta_coef_df$p_meta < 0.05) # 0
sum(meta_coef_df$q_meta < 0.05) # 0
sum(meta_coef_df$bonferroni < 0.05) # 0

write.csv(meta_coef_df, file = "Output/CpG_stress_meta_JHS_MESA_het_probes_FDR.csv")

