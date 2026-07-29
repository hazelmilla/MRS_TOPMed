# Generate EWAS catalogue for 13 LE-related CpGs at Bonferroni p < 0.05

devtools::install_github("MRCIEU/ewascatalog-r")

library(ewascatalog)
library(tidyverse)

sle_cpgs <- read.csv("~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/000 Zannas_lab/MRS TOPMed/meta_coef_sig_annotated_SLEs.csv",
                     header = TRUE) # get SLE-related CpGs

sle_cpgs <- filter(sle_cpgs, bonferroni < 0.05) # filter for 13 bonf significant sites
dim(sle_cpgs) # [1] 13 15

catalogue_results <- map_dfr(sle_cpgs$IlmnID, function(cpg){
  
  message("Querying ", cpg)
  
  out <- tryCatch(
    ewascatalog(cpg, "cpg"),
    error = function(e) NULL
  )
  
  if (is.null(out)) return(NULL)
  
  out$Query_CpG <- cpg
  
  out
}) # Get EWAS catalogue results

# convert trait names to shorten and remove redundancies

# See which names are longer than 10 characters long
idx <- which(nchar(catalogue_results$trait) > 10)
for (i in 1:length(idx)){
  num <- idx[i]
  long_names <- trimws(catalogue_results$trait[num])
  print(unique(long_names))
}
  
ewas_catalogue$traits # see all trait names

# Convert trait names
cat_results_clean <- catalogue_results %>% mutate(trait = trimws(trait)) %>%
  mutate(trait = case_when(
    trait == "C-reactive protein (CRP) levels" ~ "CRP",
    trait == "age" ~ "Age",
    trait == "Ageing" ~ "Age",
    trait == "Age 4 vs age 0" ~ "Age",
    trait == "Age" ~ "Age",
    trait == "Prevalent COPD (Self-report)" ~ "COPD",
    trait == "Incident COPD" ~ "COPD",
    trait == "Incident Lung Cancer" ~ "Lung Cancer",
    trait == "HIV infection" ~ "HIV",
    trait == "post-traumatic stress disorder" ~ "PTSD",
    trait == "Alcohol consumption per day" ~ "Alcohol",
    trait == "Tobacco smoking" ~ "Smoking",
    trait == "Prevalent Chronic Pain (Self-report)" ~ "Chronic Pain",
    trait == "Body mass index" ~ "BMI",
    trait == "Incident Type 2 Diabetes" ~ "T2 Diabetes",
    trait == "Prevalent Type 2 Diabetes (Self-report)" ~ "T2 Diabetes",
    trait == "HDL cholesterol" ~ "HDL",
    trait == "Acute myocardial infarction" ~ "MI",
    trait == "TNFRSF17 protein levels (SeqId = 2665-26)" ~ "TNFRSF17 protein",
    trait == "GZMK protein levels (SeqId = 9545-156)" ~ "GZMK protein",
    trait == "GZMA protein levels (SeqId = 3440-7)" ~ "GZMK protein",
    trait == "LTA protein levels (SeqId = 3505-6)" ~ "LTA protein",
    trait == "IL2RB protein levels (SeqId = 9343-16)" ~ "IL2RB protein",
    trait == "LY9 protein levels (SeqId = 3324-51)" ~ "LY9 protein",
    trait == "FCER2 protein levels (SeqId = 3291-30)" ~ "FCER2 protein",
    trait == "Maternal body mass index" ~ "Maternal BMI",
    trait == "maternal pre-pregnancy body mass index" ~ "Maternal BMI",
    trait == "Diet Quality Alternative Healthy Eating Index 2010 (AHEI-2010)" ~ "Diet",
    trait == "Incident Ischemic Heart Disease" ~ "Ischemic Heart Disease",
    .default = trait
    ))

# Check long trait names in cleaned dataset
idx <- which(nchar(cat_results_clean$trait) > 10)
for (i in 1:length(idx)){
  num <- idx[i]
  long_names <- trimws(cat_results_clean$trait[num])
  print(unique(long_names))
}

# Get 1 row with all traits for each CpG - make sure they are not repeated
catalogue_one_row <-
  cat_results_clean %>%
  group_by(Query_CpG) %>%
  summarise(
    n_traits = n(),
    traits = paste(unique(trait), collapse = "; "),
    .groups = "drop"
  )

# Merge with SLE EWAS data
ewas_catalogue <-
  sle_cpgs[-1] %>%
  left_join(catalogue_one_row,
            by = c("IlmnID" = "Query_CpG"))

# Write file
path = "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/000 Zannas_lab/MRS manuscript/"
write.csv(ewas_catalogue, paste0(path, "supp_table_EWAS_catalogue.csv"), row.names = FALSE)

