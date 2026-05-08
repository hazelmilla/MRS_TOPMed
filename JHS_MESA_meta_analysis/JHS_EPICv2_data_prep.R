library(tidyverse)
library(dplyr)

path = "/Users/hazelmilla/Desktop/Zannas_lab/MRS_WHI_JHS/V2"
setwd(path)

###### Data prep #######

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
