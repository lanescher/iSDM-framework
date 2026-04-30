## ---------------------------
## Objective: 
##    - Download and clean iNat data
## 
## Input:
##
## Output: 
##    - 
##
## ---------------------------

library(tidyverse)
library(SpFut.processGBIF)
library(rgbif)


# Download all of Amphibia ----------------------
# This step requires GBIF credentials. Alternatively, you can download the 
# dataset used in our analysis from: https://doi.org/10.15468/dl.38pvuq
user <- "clanescher"
pwd <- "0xjFyIlxJoKrNh"
email <- "clanescher@gmail.com"

if (file.exists("data/gbif-raw.rds")) {
  dat <- read_rds("data/gbif-raw.rds")
} else {
  dat <- download_gbif(scientificName = "Amphibia", taxonrank = "kingdom",
                       startYear = 1994, country = "US",
                       user = user, pwd = pwd, email = email)
  write_rds(dat, file = "data/gbif-raw.rds")
}
cat(dat$citation, "\n")


table(dat$dat$species)


# Isolate GPOR and DMAR ----

sp <- c("Gyrinophilus porphyriticus", "Desmognathus marmoratus")
sp.codes <- c("GPOR", "DMAR")


for (s in 1:length(sp)) {
  dat1 <- dat$dat %>%
    filter(species == sp[s])
  dat1 <- clean_gbif(dat1)
  
  
  # Get iNat data
  dat2 <- dat1 %>%
    
    filter(institutionCode == "iNaturalist") %>%
    
    dplyr::mutate(date = substr(eventDate, 1, 10),
                  date = lubridate::as_date(date),
                  year = lubridate::year(date),
                  lat = decimalLatitude,
                  lon = decimalLongitude,
                  coord.unc = coordinateUncertaintyInMeters,
                  survey.conducted = 1,
                  count = 1,
                  data.type = "PO",
                  age = "NR",
                  individual.id = NA,
                  time.to.detect = NA,
                  species = sp.codes[s]) %>%
    
    # make site.id
    dplyr::group_by(lat, lon) %>%
    dplyr::mutate(site.id = paste0(source, dplyr::cur_group_id())) %>%
    dplyr::ungroup() %>%
    
    
    # get survey.id
    dplyr::mutate(survey.id = 1:nrow(.),
                  pass.id = 1,
                  survey.pass = paste0(survey.id, "_", pass.id)) %>%
    
    
    # select cols to keep
    dplyr::select(site.id, lat, lon, stateProvince, source, coord.unc, eventDate,
                  day, month, year, survey.conducted, survey.id, pass.id,
                  survey.pass, data.type, species, age, individual.id,
                  time.to.detect, count)
  
  write.csv(dat2, paste0("data/data-ready/", sp.codes[s], "_iNat_PO.csv"))
  
  
  # Get museum data
  dat3 <- dat1 %>%
    
    filter(basisOfRecord == "PRESERVED_SPECIMEN") %>%
    
    dplyr::mutate(date = substr(eventDate, 1, 10),
                  date = lubridate::as_date(date),
                  year = lubridate::year(date),
                  lat = decimalLatitude,
                  lon = decimalLongitude,
                  coord.unc = coordinateUncertaintyInMeters,
                  survey.conducted = 1,
                  count = 1,
                  data.type = "PO",
                  age = "NR",
                  individual.id = NA,
                  time.to.detect = NA,
                  species = sp.codes[s]) %>%
    
    # make site.id
    dplyr::group_by(lat, lon) %>%
    dplyr::mutate(site.id = paste0(source, dplyr::cur_group_id())) %>%
    dplyr::ungroup() %>%
    
    
    # get survey.id
    dplyr::mutate(survey.id = 1:nrow(.),
                  pass.id = 1,
                  survey.pass = paste0(survey.id, "_", pass.id)) %>%
    
    
    # select cols to keep
    dplyr::select(site.id, lat, lon, stateProvince, source, coord.unc, eventDate,
                  day, month, year, survey.conducted, survey.id, pass.id,
                  survey.pass, data.type, species, age, individual.id,
                  time.to.detect, count)
  
  write.csv(dat3, paste0("data/data-ready/", sp.codes[s], "_museum_PO.csv"))
  
  
}
