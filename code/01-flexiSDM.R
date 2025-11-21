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

# Load libraries ----
library(tidyverse)
library(sf)
library(nimble)
library(SpFut.flexiSDM)
library(SpFut.covariates)


# Set up model variables ----
mods <- read.csv("code/model-specs.csv") %>% filter(number %in% nums.do)

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
check.covs <- mods$check.covs[a]

Bpriordist <- mods$Bprior[a]

covs.z <- c(covs.lin, covs.quad)
if ("" %in% covs.z) {
  covs.z <- covs.z[-which(covs.z == "")]  
}
if (NA %in% covs.z) {
  covs.z <- covs.z[-which(is.na(covs.z))]
}




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

data.dir <- paste0("data/", number, "_", sp.code, "_", model, "/")
if (dir.exists(data.dir) == F) {
  dir.create(data.dir)
}


# Make region ----

# Generating region requires files that are too big for github, just read it in
if (file.exists(paste0(data.dir, "region.rds"))) {
  region <- read_rds(file = paste0(data.dir, "region.rds"))
} else {
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
  
  write_rds(region, file = paste0("data/", nums.do, "_", sp.code, "_", model, "/region.rds"))
}



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
allfiles <- read.csv("data/00-data-summary-flexiSDM.csv") %>%
  filter(Species == sp.code) %>%
  rename(file.name = Data.Swamp.file.name,
         file.label = Name,
         covar.mean = Covar.mean,
         covar.sum = Covar.sum,
         data.type = Type.true) %>%
  select(file.name, file.label, covar.mean, covar.sum, data.type, PO.extent)


species.data <- load_species_data(sp.code,
                                  sp.code.all,
                                  file.info = allfiles,
                                  file.path = "data/data-ready/",
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

if (sp.code == "RACA") {
  
  if (file.exists(paste0(data.dir, "covariates.rds"))) {
    covar <- read_rds(paste0(data.dir, "covariates.rds"))
  } else {
    
    # Note that elevation data must be downloaded before running get_elevation()
    # See documentation for details.
    tri <- get_elevation(locs = region$sp.grid, path = "../species-futures/data/USA/", id.label = "conus.grid.id")
    waterbody <- get_waterbodies(locs = region$sp.grid, path = "../species-futures/data/USA/", id.label = "conus.grid.id")
    
    footprint <- get_footprint(locs = region$sp.grid, id.label = "conus.grid.id")
    climate <- get_climate(locs = region$sp.grid, id.label = "conus.grid.id")
    traveltime <- get_traveltime(locs = region$sp.grid, id.label = "conus.grid.id")
    
    covar <- full_join(tri, footprint, by = "conus.grid.id") %>%
      full_join(climate, by = "conus.grid.id") %>%
      full_join(waterbody, by = "conus.grid.id") %>%
      full_join(traveltime, by = "conus.grid.id") %>%
      mutate(sqrtarea_small = sqrt(area_small),
             sqrtarea_medium = sqrt(area_medium)) %>%
      select(conus.grid.id, sqrtarea_small, sqrtarea_medium, footprint, TRI, tmin, traveltime)
    
    write_rds(covar, file = paste0(data.dir, "covariates.rds"))
    
  }
  
}


if (sp.code == "GPOR") {
  
  if (file.exists(paste0(data.dir, "covariates.rds"))) {
    covar <- read_rds(paste0(data.dir, "covariates.rds"))
  } else {
    
    # Note that elevation and flowlines data must be downloaded before running 
    # get_elevation() and get_flowlines(). See documentation for details.
    tri <- get_elevation(locs = region$sp.grid, path = "../species-futures/data/USA/", id.label = "conus.grid.id")
    stream <- get_flowlines(locs = region$sp.grid, path = "../species-futures/data/USA/", id.label = "conus.grid.id")
    landcover <- get_landcover(locs = region$sp.grid, path = "../species-futures/data/USA/", id.label = "conus.grid.id")
    
    climate <- get_climate(locs = region$sp.grid, id.label = "conus.grid.id")
    traveltime <- get_traveltime(locs = region$sp.grid, id.label = "conus.grid.id")
    
    
    
    #### ADD N.INAT!!!!!
    
    
    
    
    covar <- full_join(tri, stream, by = "conus.grid.id") %>%
      full_join(landcover, by = "conus.grid.id") %>%
      full_join(climate, by = "conus.grid.id") %>%
      full_join(traveltime, by = "conus.grid.id") %>%
      select(conus.grid.id, streamLength.km, prec, forest, elevation, traveltime)
    
    
    covar <- covar[order(match(covar$conus.grid.id, region$sp.grid$conus.grid.id)),]
    
    write_rds(covar, file = paste0(data.dir, "covariates.rds"))
    
    
  }
}



rm <- which(complete.cases(covar[,covs.z]) == F)
if (length(rm) > 0) {
  covar <- covar[-rm,]
  region$sp.grid <- region$sp.grid[-rm,]
}

# Scale covariates 
covar_unscaled <- covar
numcols <- sapply(covar, is.numeric)
numcols <- which(numcols)
covar[,numcols] <- sapply(covar[,numcols], scale_this)



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



### Add quadratic covariates ----
if (length(covs.quad) > 0 & paste0(covs.quad, collapse = "") != "") {
  for (c in 1:length(covs.quad)) {
    covar[,paste0(covs.quad[c], "2")] <- covar[,covs.quad[c]] * covar[,covs.quad[c]]
    covs.z <- c(covs.z, paste0(covs.quad[c], "2"))
  }
}



# NIMBLE ----

sp.data <- sppdata_for_nimble(species.data,
                              region,
                              file.info = allfiles,
                              covar = covar,
                              covs.inat = covs.inat,
                              covs.PO = covs.PO,
                              DND.maybe = 1,
                              keep.conus.grid.id = gridkey$conus.grid.id[which(gridkey$group == "train")]) # get rid of PO cells that are in the wrong grid cells


### Data/constants ----
tmp <- data_for_nimble(sp.data, covar = covar, covs.z,
                       sp.auto = sp.auto, coarse.grid = coarse.grid, region = region,
                       process.intercept = F,
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
                    Bprior = Bpriordist,
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


source("../species-futures/code/03-species-models/xx-flexiSDM-setuptests.R")


# Remove local and block in case the setup is run locally but the model is fit on the HPC.
# Remove other unnecessary files to reduce the size of setup_BLOCK.Rdata
rm(list=c('local','block','args','conus.covar.grid','conus.grid','usa','conus',
          'pl','centroid'))


# Save environment and full set up
save.image(paste0(out.dir, "setup_",block.out,".Rdata"))


# End script - proceed to 02-flexiSDM.R
