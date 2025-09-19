## ---------------------------
## Objective: 
##    - Imports species and covariate data
##    - Scales covariates
##    - Sets up data for NIMBLE
## 
## Input:
##    - functions/FXN-MVPv1.1.R
##    - code/03-species-models/MVPv1.csv
##
## Output: 
##    - setup_BLOCK.rdata (saves environment for import to fit model)
##
## ---------------------------

start1 <- Sys.time() 
print(paste0('Beginning 01-flexiSDM script at ', start1))


# summarize datasets to get up-to-date available datasets
# this script also prints messages from QC functions
# source("DATA SWAMP/00-summarize-datasets.R")


# EDIT THIS SECTION ----
nums.do <- 4
block <- "none"
# block <- c("none", 1, 2, 3)
local <- 1
a <- 1
# ---


# DO NOT EDIT BELOW HERE ----
args = commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  
  # Number of model to run
  nums.do = as.numeric(args[1])
  
  # Assign a block to run if running locally
  block = as.numeric(args[2])
  
  if (block == 4) {
    block <- 'none'
  }  
  
  # Running locally? Yes = 1
  local = as.numeric(args[3])
  
  if (local == 0) {
    setwd('/caldera/hovenweep/projects/usgs/ecosystems/eesc/cscher/iSDM-framework/')
  } 
} 

# load libraries ----
library(tidyverse)
library(sf)
library(nimble)
library(SpFut.flexiSDM)



# Set up model variables ----
mods <- read.csv("code/MVPv1.csv") %>% filter(number %in% nums.do)

tmp <- expand.grid(block.out = block, number = unique(mods$number))
mods <- full_join(mods, tmp, by = c("number"))



# Get variables for model from MVPv1.csv
number <- mods$number[a]
sp.code <- mods$sp.code[a]
model <- mods$model[a]

sp.auto <- mods$sp.auto[a]
coarse.grid <- mods$coarse.grid[a]
year.start <- mods$year.start[a]
year.end <- mods$year.end[a]
buffer <- mods$buffer[a]
filter.region <- mods$filter.region[a]
spat.bal <- mods$spat.bal[a]
coordunc <- mods$coordunc[a]
coordunc_na.rm <- mods$coordunc_na.rm[a]
block.folds <- mods$block.folds[a]
block.rows <- mods$block.rows[a]
block.cols <- mods$block.cols[a]
block.out <- mods$block.out[a]

if (block.out == "none") {
  blockname <- "full"
} else {blockname <- block.out}

zero_mean <- mods$zero_mean[a]
tau <- mods$tau[a]


iter <- mods$iter[a]
thin <- mods$thin[a]
burnin <- floor(iter*0.75)
region.sub <- mods$region.sub[a]
lat.hi <- mods$lat.hi[a]
lat.lo <- mods$lat.lo[a]
lon.hi <- mods$lon.hi[a]
lon.lo <- mods$lon.lo[a]

project <- mods$project[a]
if (block.out != "none") {
  project <- 0
}

# Get covariates from MVPv1.csv
covs.PO <- unlist(str_split(mods$covs.PO[a], pattern = ", "))
covs.inat <- unlist(str_split(mods$covs.inat[a], pattern = ", "))
covs.lin <- unlist(str_split(mods$covs.lin[a], pattern = ", "))
covs.quad <- unlist(str_split(mods$covs.quad[a], pattern = ", "))
covs.int.factor <- unlist(str_split(mods$covs.int.factor[a], pattern = ", "))
reference <- mods$reference[a]
covs.int.cont <- unlist(str_split(mods$covs.int.cont[a], pattern = ", "))
check.covs <- mods$check.covs[a]

Bpriordist <- mods$Bpriordist[a]
Bpriorvar1 <- mods$Bpriorvar1[a]
Bpriorvar2 <- mods$Bpriorvar2[a]

process.intercept <- F


codeKey <- read.csv("data/model-specieslist.csv")
common <- codeKey %>%
  filter(DS.code == sp.code) %>%
  pull(SSAR.common)
sp.code.all <- codeKey %>%
  filter(DS.code == sp.code) %>%
  pull(all.codes)
sciname <- codeKey %>%
  filter(DS.code == sp.code) %>%
  pull(SSAR.scientific)

gen <- stringr::word(sciname, 1)

if (gen %in% c("Acris", "Anaxyrus", "Aquarana", "Ascaphus", 
               "Craugastor", "Dendrobates", "Dryophytes", "Eleutherodactylus",
               "Gastrophryne", "Glandirana", "Hyla", "Hypopachus", "Incilius",
               "Leptodactylus", "Lithobates", "Osteopilus",
               "Pseudacris", "Rana", "Rhinella", "Rhinophrynus",
               "Scaphiopus", "Smilisca", "Spea", "Xenopus")) spp.type <- "Frog/Toad"
if (gen %in% c("Ambystoma", "Amphiuma", "Aneides", "Batrachoseps",
               "Cryptobranchus", "Desmognathus", "Dicamptodon", 
               "Ensatina", "Eurycea", "Gyrinophilus", "Hemidactylium",
               "Hydromantes", "Necturus", "Notophthalmus",
               "Phaeognathus", "Plethodon", "Pseudobranchus",
               "Pseudotriton", "Rhyacotriton", "Siren", 
               "Stereochilus", "Taricha", "Urspelerpes")) spp.type <- "Salamander"




# Ready to start processing ----
print(paste0("MODEL: ",number,"_",sp.code,"_",model))
print(paste0("BLOCK: ", block.out))
print(paste0("PROJECTIONS: ", project))
print(paste0("COVARIATES: ", paste0(c(covs.lin, covs.quad), collapse = ", ")))


# Make output folder ----
out.dir <- paste0("outputs/", number, "_", sp.code, "_", model, "/")
if (dir.exists(out.dir) == F) {
  dir.create(out.dir)
}


# Make region ----
# Load ranges
range.path <- c(paste0("data/", sp.code, "/GAP/"),
                paste0("data/", sp.code, "/IUCN/"))
range.name <- c("GAP", "IUCN")
rangelist <- get_range(range.path,
                       range.name,
                       crs = 4326)


# USA boundary
exclude <- c("Alaska", "Hawaii", "Commonwealth of the Northern Mariana Islands",
             "American Samoa", "United States Virgin Islands", "Guam", "Puerto Rico")
usa <- st_read("../species-futures/data/USA/maps/cb_2018_us_state_500k/cb_2018_us_state_500k.shp") %>%
        filter((NAME %in% exclude) == F) %>%
        st_union()
# CONUS grid
load("../species-futures/data/USA/grid-covar.rdata")

# Region
region <- make_region(rangelist,
                      buffer = buffer,
                      crs = 3857,
                      sub = region.sub,
                      boundary = usa,
                      grid = conus.grid,
                      rm.clumps = T,
                      clump.size = 50)

# Set up cross validation blocks ----
spatblocks <- make_CV_blocks(region, rows = block.rows, cols = block.cols, k = block.folds)




# Set up training and test sets based on cross validation blocks
if (block.out == "none") {
  test.i <- c()
  train.i <- region$sp.grid$conus.grid.id
  
} else {
  block1 <- spatblocks %>%
    filter(folds == block.out)
  
  # find grid.ids for test block, everything else is train
  test.i1 <- st_intersection(region$sp.grid, block1)
  test.i <- st_intersection(region$sp.grid, block1) %>%
    pull(conus.grid.id) %>%
    unique()
  train.i <- filter(region$sp.grid, conus.grid.id %in% test.i == F) %>%
    pull(conus.grid.id)
  
}



# Make gridkey and spatkey ----
gridkey <- select(region$sp.grid, conus.grid.id) %>%
  st_drop_geometry() %>%
  mutate(grid.id = 1:nrow(.),
         group = case_when(conus.grid.id %in% train.i ~ "train",
                           conus.grid.id %in% test.i ~ "test"))

if (coarse.grid == T) {
  spatRegion <- suppressWarnings(make_spatkey(region$sp.grid))
  
  # Print spatial grid
  if (block.out == 'none') {
    pl <- ggplot(spatRegion$spat.grid) + geom_sf() + theme_bw()
    ggsave(pl, file = paste0(out.dir, "2_inputmap-e_spatGrid.jpg"), height = 8, width = 10)
  }
} else {
  spatRegion <- NULL
}



# Species data ----
# get all files that have data for that species
allfiles <- read.csv("data/dataset-summary-full.csv")
tmp1 <- gsub("[|]", "$|^", sp.code.all)
tmp2 <- paste0("^", tmp1, "$")
allfiles <- allfiles[grep(tmp2, allfiles$species),] %>%
  select(-species, -percentdet) %>%
  distinct()

# detection covariates for each of these datasets
covs <- read.csv("data/00-data-summary-flexiSDM.csv") %>%
  filter(Data.Swamp.file.name %in% allfiles$file) %>%
  select(Data.Swamp.file.name, Covar.mean, Covar.sum)
covs <- covs[order(match(covs$Data.Swamp.file.name, allfiles$file)),]

covariates <- list()
for (i in 1:nrow(covs)) {
  covs.mean <- unlist(strsplit(covs$Covar.mean[i], split = ", "))
  covs.sum <- unlist(strsplit(covs$Covar.sum[i], split = ", "))
  #area <- unlist(strsplit(covs$Area[i], split = ","))
  covs1 <- c(covs.mean, covs.sum)
  covs1 <- covs1[which(is.na(covs1) == F)]
  covariates[[covs$Data.Swamp.file.name[i]]] <- covs1
}

species.data <- load_species_data(sp.code,
                                  sp.code.all,
                                  file.name = allfiles$file,
                                  file.label = allfiles$name,
                                  file.path = "data/data-ready/",
                                  keep.cols = covariates,
                                  region = region, 
                                  filter.region = filter.region,
                                  year.start = year.start,
                                  year.end = year.end,
                                  coordunc = coordunc,
                                  coordunc_na.rm = coordunc_na.rm,
                                  spat.thin = spat.bal,
                                  keep.conus.grid.id = gridkey$conus.grid.id[which(gridkey$group == "train")])





### Plot species data ----
if (block == "none") {
  title <- ", full model"
} else {
  title <- paste0(", excluding block ", block)
}

pl <- map_species_data(region = region,
                       species.data = species.data,
                       year.start = year.start,
                       year.end = year.end,
                       plot = "samples",
                       plot.region = T,
                       details = T,
                       title = paste0(common, " (", sp.code, ")", title))
ggsave(pl, file = paste0(out.dir, "2_inputmap-b_data-", blockname, "-details.jpg"), height = 8, width = 10)

pl <- map_species_data(region = region,
                       species.data = species.data,
                       year.start = year.start,
                       year.end = year.end,
                       plot = "samples",
                       plot.region = T,
                       details = F,
                       title = paste0(common, " (", sp.code, ")", title))
ggsave(pl, file = paste0(out.dir, "2_inputmap-a_data-", blockname, ".jpg"), height = 8, width = 10)


if (block.out != "none") {
  
  pl <- map_species_data(region = region,
                         species.data = species.data,
                         year.start = year.start,
                         year.end = year.end,
                         plot = "samples",
                         plot.blocks = T,
                         blocks = spatblocks[which(spatblocks$folds == block),],
                         plot.region = T,
                         details = T,
                         title = paste0(common, " (", sp.code, ")", title))
  ggsave(pl, file = paste0(out.dir, "2_inputmap-d_blocks-", blockname, "-details.jpg"), height = 8, width = 10)
  
  pl <- map_species_data(region = region,
                         species.data = species.data,
                         year.start = year.start,
                         year.end = year.end,
                         plot = "samples",
                         plot.blocks = T,
                         blocks = spatblocks[which(spatblocks$folds == block),],
                         plot.region = T,
                         details = F,
                         title = paste0(common, " (", sp.code, ")", title))
  ggsave(pl, file = paste0(out.dir, "2_inputmap-c_blocks-", blockname, ".jpg"), height = 8, width = 10)
  
  
} else {
  pl <- map_species_data(region = region,
                         species.data = species.data,
                         year.start = year.start,
                         year.end = year.end,
                         plot = "samples",
                         plot.blocks = T,
                         blocks = spatblocks,
                         plot.region = T,
                         details = T,
                         title = paste0(common, " (", sp.code, ")", title))
  ggsave(pl, file = paste0(out.dir, "2_inputmap-d_blocks-", blockname, "-details.jpg"), height = 8, width = 10)
  
  pl <- map_species_data(region = region,
                         species.data = species.data,
                         year.start = year.start,
                         year.end = year.end,
                         plot = "samples",
                         plot.blocks = T,
                         blocks = spatblocks,
                         plot.region = T,
                         details = F,
                         title = paste0(common, " (", sp.code, ")", title))
  ggsave(pl, file = paste0(out.dir, "2_inputmap-c_blocks-", blockname, ".jpg"), height = 8, width = 10)
  
}





# Covariate data ----
load("../species-futures/data/USA/grid-covar.rdata")
covar <- conus.covar.grid %>%
  filter(conus.grid.id %in% region$sp.grid$conus.grid.id)
covar <- covar[order(match(covar$conus.grid.id, region$sp.grid$conus.grid.id)),]

# load iNat data for similar species
if ("iNaturalist" %in% names(species.data$obs)) {
  cat("loading iNat records of supplemental species\n")
  
  
  
  # Get basic covariate layer
  inat <- readRDS("data/USA/inateffortcov.rds") %>%
    mutate(conus.grid.id = as.character(conus.grid.id))
  
  covar <- dplyr::left_join(covar, inat, by = c("conus.grid.id"))
  
  # Assign Frog or Salamander iNat counts
  if (spp.type == 'Frog/Toad') {
    covar$n.inat <- covar$n.frog
  } else {
    covar$n.inat <- covar$n.sal
  }
  covar <- select(covar, -n.frog, -n.sal)
  
  
  
  # Now check if records for this species need to be added
  
  # read in all inat data for this species
  files <- list.files("../species-futures/DATA SWAMP/data-ready/", pattern = "_iNat_PO.csv")
  fs <- files[grep(sp.code.all, files)]
  if (length(fs) == 0) next
  inat.dat1 <- c()
  for (f in 1:length(fs)) {
    
    dat <- read.csv(paste0("../species-futures/DATA SWAMP/data-ready/", fs[f])) %>%
      filter(year >= year.start, year <= year.end)
    
    # iNaturalist records for subspecies names might be duplicates, remove here
    new <- paste0(dat$lat, dat$lon, dat$day, dat$month, dat$year)
    old <- paste0(inat.dat1$lat, inat.dat1$lon, inat.dat1$day, inat.dat1$month, inat.dat1$year)
    
    if (length(old) > 0) {
      rm <- which(new %in% old)
      dat <- dat[-rm,]
      
      if (nrow(dat) == 0) {
        #cat("All iNaturalist records are potential duplicates\n")
        next
      } else {
        inat.dat1 <- bind_rows(inat.dat1, dat)
      }
    } else {
      inat.dat1 <- bind_rows(inat.dat1, dat)
    }
  } # end loading inat files
  
  
  # if all points in a state are obscured (coord.unc > 25000), species was not included
  tmp <- inat.dat1 %>%
    group_by(stateProvince) %>%
    summarize(min.unc = min(coord.unc, na.rm = T)) %>%
    mutate(min.unc = case_when(min.unc == Inf ~ NA,
                               T ~ min.unc)) %>%
    filter(min.unc > 25000)
  obsc.state <- tmp$stateProvince
  statecodes <- read.csv("data/statecodes.csv")
  obsc.state <- statecodes$abbrev[which(statecodes$state %in% obsc.state)]
  
  # species was already included
  if (nrow(tmp) == 0) {
    
    
    # Need to add species
  } else {
    cat("Adding", sp.code, "to iNat effort covariate\n")
    
    inat.cells <- inat.dat1 %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      st_transform(crs = 3857) %>%
      st_intersection(region$sp.grid) %>%
      st_drop_geometry() %>%
      group_by(conus.grid.id) %>%
      summarize(n.inat.sp = n()) %>%
      full_join(region$sp.grid, ., by = "conus.grid.id") %>%
      st_drop_geometry()
    inat.cells$n.inat.sp[is.na(inat.cells$n.inat.sp)] <- 0
    
    covar <- full_join(covar, select(inat.cells, conus.grid.id, n.inat.sp), by = 'conus.grid.id') %>%
      mutate(n.inat = n.inat + n.inat.sp)
  }
}





# add centroid lat and lon of each grid cell to covar
centroid <- st_centroid(region$sp.grid) %>%
              st_coordinates() %>%
              as.data.frame() %>%
              mutate(conus.grid.id = region$sp.grid$conus.grid.id)

colnames(centroid)[1:2] <- c("lon", "lat")
covar <- left_join(covar, centroid, by = "conus.grid.id")



### Read in process covariates ----
covs.z <- c(covs.lin, covs.quad, covs.int.factor)
if ("" %in% covs.z) {
  covs.z <- covs.z[-which(covs.z == "")]  
}
if (NA %in% covs.z) {
  covs.z <- covs.z[-which(is.na(covs.z))]
}



### Add additional/derived covariates -----
# If you add a new covariate, add a row to data/covariate-labels.csv with a label
if (sp.code == "RACA") {
  
  covar <- covar %>%
    mutate(sqrtarea_small = sqrt(area_small),
           sqrtarea_medium = sqrt(area_medium))
  
} else if (sp.code == "PLSE") {
  
  covar <- covar %>%
    mutate(N_all = N + NW + NE)
  
}

rm <- which(complete.cases(covar[,covs.z]) == F)
if (length(rm) > 0) {
  covar <- covar[-rm,]
  region$sp.grid <- region$sp.grid[-rm,]
}

### Scale covariates ----
covar_unscaled <- covar

# scale numeric cols
numcols <- sapply(covar, is.numeric)
numcols <- which(numcols)
covar[,numcols] <- sapply(covar[,numcols], scale_this)

# remove covariates that are correlated > 0.4
if (check.covs == T){
  covs.rm <- select_covar(covs.z, threshold = 0.4)
  covs.lin <- covs.lin[-which(covs.lin %in% covs.rm)]
  covs.quad <- covs.quad[-which(covs.quad %in% covs.rm)]
  
  covs.z <- c(covs.lin, covs.quad)
}

if (length(covs.z) < 3) {
  stop("There are fewer than 3 covariates remaining! This probably isn't a good model")
}


### Plot covariates ----
if (block.out == "none") {
  
  
  # Process covs
  
  # get covariate labels
  covlabs <- read.csv("data/covariate-labels.csv") %>%
    filter(covariate %in% covs.z)
  
  plot_covar(covar,
             region,
             cov.names = covlabs$covariate,
             cov.labels = covlabs$Label,
             out.path = out.dir,
             out.name = "1_covariates-a_process-map")
  
  cor_covar(covar, 
            cov.names = covlabs$covariate,
            cov.labels = covlabs$Label,
            out.path = out.dir,
            out.name = "1_covariates-a_process-correlations", 
            color.threshold = 0.25)
  
  
  
  # iNat covs
  if ("iNaturalist" %in% names(species.data$obs)) {
    # get covariate labels
    covlabs <- read.csv("data/covariate-labels.csv") %>%
      filter(covariate %in% covs.inat)
    
    plot_covar(covar,
               region,
               cov.names = covlabs$covariate,
               cov.labels = covlabs$Label,
               out.path = out.dir,
               out.name = "1_covariates-b_iNat-map")
    
    if (length(covs.inat) > 1) {
      cor_covar(covar, 
                cov.names = covlabs$covariate,
                cov.labels = covlabs$Label,
                out.path = out.dir,
                out.name = "1_covariates-b_iNat-correlations", 
                color.threshold = 0.25)
    }
  }
  
  
  
  
  # PO covs
  
  # get covariate labels
  covlabs <- read.csv("data/covariate-labels.csv") %>%
    filter(covariate %in% covs.PO)
  
  plot_covar(covar,
             region,
             cov.names = covlabs$covariate,
             cov.labels = covlabs$Label,
             out.path = out.dir,
             out.name = "1_covariates-c_PO-map")
  
  if (length(covs.PO) > 1) {
    cor_covar(covar, 
              cov.names = covlabs$covariate,
              cov.labels = covlabs$Label,
              out.path = out.dir,
              out.name = "1_covariates-c_PO-correlations", 
              color.threshold = 0.25)
  }
  
}



### Add quadratic covariates and interactions ----
if (length(covs.quad) > 0 & paste0(covs.quad, collapse = "") != "") {
  for (c in 1:length(covs.quad)) {
    covar[,paste0(covs.quad[c], "2")] <- covar[,covs.quad[c]] * covar[,covs.quad[c]]
    covs.z <- c(covs.z, paste0(covs.quad[c], "2"))
  }
}

# and interactions
if (length(covs.int.factor) == 1 & is.na(covs.int.factor) == F) {
  
  # pivot factor wider to get dummy cols
  covar[,covs.int.factor] <- paste0(covs.int.factor, ".", covar[,covs.int.factor])
  covar <- covar %>%
    mutate(fact = 1) %>%
    pivot_wider(names_from = all_of(covs.int.factor), values_from = fact, values_fill = 0)
  
  covs.z <- c(covs.z, colnames(covar)[grep(paste0(covs.int.factor, "."), colnames(covar))])
  
  for (j in 1:length(covs.int.cont)) {
    # add interaction
    tmp <- add_int_cols(covar, int.factor = covs.int.factor, int.cont = covs.int.cont[j])
    covar <- tmp$covar
    covs.z <- c(covs.z, tmp$names)
  }
  
  # remove reference level
  rm <- grep(reference, covs.z)
  covs.z <- covs.z[-rm]
  
  # remove original column name
  rm <- grep(covs.int.factor, covs.z)[1]
  covs.z <- covs.z[-rm]
}











# NIMBLE ----

summary <- read.csv("data/00-data-summary-flexiSDM.csv") %>%
  filter(Data.Swamp.file.name %in% allfiles$file,
         Name %in% names(species.data$obs)) %>%
  select(-Data.Swamp.file.name) %>%
  distinct()
summary <- summary[order(match(summary$Name, names(species.data$obs))),] %>% distinct()
summary$Area <- "" # we're not using the area column anymore


sp.data <- sppdata_for_nimble(species.data,
                              region,
                              data.type = summary$Type.true,
                              PO.extent = summary$PO.extent,
                              covar = covar,
                              covs.inat = covs.inat,
                              covs.PO = covs.PO,
                              covs.mean = summary$Covar.mean,
                              covs.sum = summary$Covar.sum,
                              offset.area = summary$Area,
                              DND.maybe = 1,
                              keep.conus.grid.id = gridkey$conus.grid.id[which(gridkey$group == "train")]) # get rid of PO cells that are in the wrong grid cells


### Data/constants ----
tmp <- data_for_nimble(sp.data, covar = covar, covs.z,
                       sp.auto = sp.auto, coarse.grid = coarse.grid, region = region,
                       process.intercept = process.intercept,
                       gridkey = gridkey, spatRegion= spatRegion)

data <- tmp$data
constants <- tmp$constants


# Add state indicator variable for iNat data to indicate which states have taxon geoprivacy
# Add state indicator for multi-state PO to indicate which states have data
constants <- add_state_ind(species.data,
                           region,
                           gridkey,
                           constants,
                           covs.inat = covs.inat,
                           obsc.state = obsc.state,
                           keep.conus.grid.id = gridkey$conus.grid.id[which(gridkey$group == "train")])




### Code ----
code <- nimble_code(data,
                    constants, 
                    path = out.dir,
                    sp.auto = sp.auto, 
                    coarse.grid = coarse.grid,
                    Bpriordist = Bpriordist, Bpriorvar1 = Bpriorvar1, Bpriorvar2 = Bpriorvar2,
                    block.out = block.out,
                    min.visits.incl = 3, 
                    zero_mean = zero_mean,
                    rm.state = F,
                    tau = tau)

### Initial values ----
inits <- function(x){nimble_inits(data,
                                  constants,
                                  sp.auto = sp.auto,
                                  seed = x)}

### Parameters ----
params <- nimble_params(data,
                        constants,
                        lambda = T,
                        XB = T,
                        sp.auto = sp.auto,
                        effort = T)



end1 <- Sys.time() - start1

# Remove local and block in case the setup is run locally but the model is fit on the HPC.
# Remove other unnecessary files to reduce the size of setup_BLOCK.Rdata
rm(list=c('local','block','args','conus.covar.grid','conus.grid','usa','conus',
          'pl','centroid'))


# Save environment and full set up
save.image(paste0(out.dir, "setup_",block.out,".Rdata"))


# End script - proceed to 02-flexiSDM.R
