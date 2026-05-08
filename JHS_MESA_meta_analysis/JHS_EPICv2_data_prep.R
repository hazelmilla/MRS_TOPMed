##### Data prep for MRS analysis in EPICv2 Data #####
library(dplyr)
library(tidyverse)

path = "/work/users/h/e/hemilla/WHI_JHS_MRS/EPICv2"
setwd(path)

data <- read.csv("analysis1_2020_update_dbgap_id.csv", header = TRUE)
pc <- read.delim("JHS_pcair_pcs_output_TOEID.txt", header = TRUE, sep = "\t")
horvath <- read.csv("DNAmAgeCalcProject_14998_Results.csv", header = TRUE)
marital <- read.csv("marriage_status_dbgap_id.csv")
subjid <- read.csv("subjid_TOE.csv", header = TRUE)
mrs <- read.csv("MRS_JHS_EPICV2.csv", header = TRUE)

# get pc & horvath outputs together for genetic ancestry pcs 
# and Houseman cell type proportions
colnames(horvath)[1] = "TOEID"
pc_horvath <- merge(pc, horvath, by = "TOEID")
pc_ctp <- merge(pc_horvath, subjid, by = "TOEID")

# merge with data and marital to have all phenotypes in 1 place
pheno <- merge(data, marital, by = "SUBJECT_ID")
dat <- merge(pheno, pc_ctp, by = "SUBJECT_ID")

# subset the data 
dat <- dat %>% subset(select = c("SUBJECT_ID", "perceivedStress", "sex", "age",
                                 "VisitDate", "alc", "currentSmoker", "BMI", 
                                 "Income", "edu3cat", "nbSESanascore",
                                 "SOCA1A", "PC1", "PC2", "PC3", "PC4", "PC5", 
                                 "PC6", "PC7", "PC8", "PC9", "PC10", "NK", 
                                 "Bcell", "Mono", "Gran", "CD8T", "CD4T"))

# merge with MRS scores
data <- merge(dat, mrs, by = "SUBJECT_ID")

med <- median(data$perceivedStress, na.rm = TRUE)

data$stress2cat[data$perceivedStress <= med] <- 1
data$stress2cat[data$perceivedStress > med] <- 2

med_mrs <- median(data$MRS, na.rm = TRUE)

data$MRS2cat[data$MRS <= med_mrs] <- 1
data$MRS2cat[data$MRS > med_mrs] <- 2

med_mrs_FDR <- median(data$MRS_FDR, na.rm = TRUE)

data$MRS2cat_FDR[data$MRS_FDR <= med_mrs_FDR] <- 1
data$MRS2cat_FDR[data$MRS_FDR > med_mrs_FDR] <- 2

med_mrs_bonf <- median(data$MRS_bonferroni, na.rm = TRUE)

data$MRS2cat_bonf[data$MRS_bonferroni <= med_mrs_bonf] <- 1
data$MRS2cat_bonf[data$MRS_bonferroni > med_mrs_bonf] <- 2

data$marital2cat[data$SOCA1A == "N"] <- "Not_married"
data$marital2cat[data$SOCA1A == "M"] <- "Married"
data$marital2cat[data$SOCA1A == "D"] <- "Not_married"
data$marital2cat[data$SOCA1A == "S"] <- "Not_married"
data$marital2cat[data$SOCA1A == "W"] <- "Not_married"
table(data$marital2cat)

# save file
write.csv(data, "prepared_data_v2.csv", row.names = FALSE)

###### CHD data prep #######
path = "/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/V2"
setwd(path)

data_unfiltered <- read.csv("prepared_data_v2.csv", header = TRUE)
dim(data_unfiltered)
#[1] 1617   36

analysis <- read.csv("analysis1_2020_update_dbgap_id.csv", 
                     header = TRUE)

chd <- read.csv("incevtchd_2020_update_dbgap_id.csv", header = TRUE)
dim(chd)
#[1] 3406   14

# Data prep
chd_status <- analysis %>% subset(select = c("SUBJECT_ID", "CHDHx"))
chd_outcome <- chd %>% subset(select = c("SUBJECT_ID", "CHD", "V1date", "date"))

CHDHx <- filter(chd_status, CHDHx == 1) # no of participants = 256
noCHDHx <- filter(chd_status, CHDHx != 1) # no of participants = 3150

# Remove participants with prior history of CHD
data1 <- data_unfiltered %>% filter(data_unfiltered$SUBJECT_ID %in% noCHDHx$SUBJECT_ID)
dim(data1)
#[1] 1499   36

#convert date formats
chd_outcome <- chd_outcome %>%
  mutate(
    date_format = as.Date(date, format = "%m/%d/%Y"),
    v1date_format = as.Date(V1date, format = "%m/%d/%Y")
  )

# Calculate the time difference in days
chd_outcome <- chd_outcome %>%
  mutate(
    time_to_event = difftime(date_format, v1date_format, units = "days")
  ) %>%
  mutate(
    days_to_event = as.numeric(time_to_event))

data_unfiltered <- merge(data1, chd_outcome)
write.csv(data_unfiltered, "chd_data_v2.csv", row.names = FALSE)
