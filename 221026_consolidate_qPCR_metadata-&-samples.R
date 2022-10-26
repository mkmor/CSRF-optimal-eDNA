######################################
##### Date Created  ## 17 Oct 22 #####
##### Date Modified ## 25 Oct 22 #####
######################################

library(readxl)
library(tidyr)
library(tidyverse)
library(dplyr)
library(lubridate) # to fix inconsistent date/time formats imported from excel



qPCRfiles=list.files(path = "/Users/morrisonme/Documents/Projects", 
                     pattern = "opt-eDNA-[0-9]{2}_qPCR_data", 
                     full.names = TRUE)

metadata = lapply(qPCRfiles, function(x) read_excel(x, sheet = 2))
samples = lapply(qPCRfiles, function(x) read_excel(x, sheet = 3))


names(metadata) = qPCRfiles
names(samples) = qPCRfiles

# Ensures dates are in unified format across datasets.
for (i in 1:length(metadata)){
  if ((typeof(metadata[[i]]$eventDate) == "double") == TRUE ) {
    metadata[[i]]$newDate = as.Date(as.character(metadata[[i]]$eventDate))
  } else {
    ifelse( 
      ((nchar(metadata[[i]]$eventDate) == 5) == TRUE),
      (metadata[[i]]$newDate = as.Date(as.numeric(metadata[[i]]$eventDate), origin = "1899-12-30")),
      (metadata[[i]]$newDate2 = as.Date(parse_date_time(metadata[[i]]$eventDate, orders = "d m y"))))
  }
  metadata[[i]]$eventDate = metadata[[i]]$newDate  # merge Date columns
  metadata[[i]]$eventDate[!is.na(metadata[[i]]$newDate2)] = metadata[[i]]$newDate2[!is.na(metadata[[i]]$newDate2)]
}

# match event date to samples. 
# Always have the vector of smaller size as second arg in match()
for (j in 1:length(samples)) {
  samples[[j]]$date = metadata[[j]]$eventDate[match(samples[[j]]$eventID, metadata[[j]]$eventID)]
  samples[[j]]$year = year(samples[[j]]$date)
  samples[[j]] = samples[[j]] %>% 
    drop_na(year) %>%
    mutate(det_prob = case_when(
      occurrenceStatus == "Detected" ~ 1.00,
      occurrenceStatus == "Suspected" ~ 0.66,
      occurrenceStatus == "Inconclusive" ~ 0.33,
      occurrenceStatus == "Not detected" ~ 0.00)) %>%
    subset(select=c(eventID, scientificName, occurrenceStatus, det_prob, year, date))
  } 

# combine all sample df into a single df for plotting
samples.df = samples %>%
  map_df(~bind_rows(.x), .id = "projectID")
