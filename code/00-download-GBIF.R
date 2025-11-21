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




# Get iNat data ----

if (file.exists("data/gbif-raw.rds")) {

  cat("Loading raw data\n")
  dat <- read_rds("data/gbif-raw.rds")
} else {

  cat("Downloading raw data\n")
  
  x <- rgbif::occ_download(
    
    rgbif::pred("hasGeospatialIssue", FALSE),
    rgbif::pred("hasCoordinate", TRUE),
    rgbif::pred("occurrenceStatus", "PRESENT"),
    rgbif::pred("taxonKey", 131),
    rgbif::pred("country", "US"),
    rgbif::pred_gte("year", 2010),
    format = "SIMPLE_CSV",
    user = "rmummah", pwd = "M@gicPizza1", email = "rileymummah@gmail.com")
  
  # get download status
  dlKey <- rgbif::occ_download_wait(x)
  
  # Load data
  dat.gbif <- rgbif::occ_download_get(dlKey$key, overwrite = T) %>% rgbif::occ_download_import()
  
  # get citation
  cite <- rgbif::gbif_citation(as.character(dlKey$key))[[1]]
  
  dat <- list(dat = dat.gbif,
              citation = cite)
  
  write_rds(dat, "data/gbif-raw.rds")
}
cat(dat$citation, "\n")


table(dat$dat$species)


# Isolate GPOR and RACA ----

sp <- c("Gyrinophilus porphyriticus", "Rana cascadae")
sp.codes <- c("GPOR", "RACA")


for (s in 1:length(sp)) {
  dat <- dat$dat %>%
    filter(species == sp[s])
  dat <- clean_gbif(dat)
  
  
  dat1 <- dat %>%
    
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
  
  
  write.csv(dat1, paste0("data/data-ready/", sp.codes[s], "_iNat_PO.csv"))
  
}
