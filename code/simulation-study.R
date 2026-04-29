

# Simulation study for model framework paper


# remotes::install_github("rileymummah/SpFut.flexiSDM-dev@v1.1-updates", force = T, auth_token = Sys.getenv("pat"))
library(SpFut.flexiSDM)
library(sf)
library(nimble)
library(tidyverse)



parall <- c()
for (s in 1) {
  
  set.seed(s)
  
  # Simulate species ----
  sp <- sim_species(grid.dim = 20, 
                    shape = "hexagon",
                    cellsize.sim = 1,         
                    cellsize.inf = 1,        
                    covForm = c("linear", "linear", "quadratic"),
                    covCorGroup = c("A", "B", "C"),
                    covType = c("spatial", "spatial", "spatial"),
                    covNoise = c(NA, NA, NA),
                    paramShift = 0,
                    paramSD = 0.1,
                    spatial = T)
  
  # Simulate data ----
  
  # PO data
  po_dat <- sim_POdat(true = sp, 
                      shape = "hexagon",
                      link = "log",
                      covForm = c("linear"),
                      covCorGroup = c("A"),
                      covType = c("spatial"),
                      covNoise = c(NA),
                      paramSD = 0.5)
  
  # count 1 visit
  count_dat1 <- sim_surveydat(true = sp,
                             survey.type = "count",
                             method = "random",
                             nsite = 10000,
                             nvisit = 1,
                             nvisit.sd = 0,
                             link = "log",
                             min.visits.incl = Inf,
                             covFormSite = c("linear"),
                             covCorGroupSite = c("A"),
                             covFormVisit = c(""),
                             covCorGroupVisit = c(""))
  
  # count 100 visits
  count_dat2 <- sim_surveydat(true = sp,
                             survey.type = "count",
                             method = "random",
                             nsite = 10,
                             nvisit = 100,
                             nvisit.sd = 0,
                             link = "logit",
                             min.visits.incl = 10,
                             covFormSite = c(""),
                             covCorGroupSite = c(""),
                             covFormVisit = c("linear"),
                             covCorGroupVisit = c("A"))
  
  # DND 1 visit
  dnd_dat1 <- sim_surveydat(true = sp,
                           survey.type = "DND",
                           method = "regular",
                           nsite = 10000,
                           nvisit = 1,
                           nvisit.sd = 0,
                           link = "log",
                           min.visits.incl = Inf,
                           covFormSite = c("linear"),
                           covCorGroupSite = c("A"),
                           covFormVisit = c(""),
                           covCorGroupVisit = c(""))
  
  
  # DND 100 visits
  dnd_dat2 <- sim_surveydat(true = sp,
                           survey.type = "DND",
                           method = "regular",
                           nsite = 10,
                           nvisit = 100,
                           nvisit.sd = 0,
                           link = "logit",
                           min.visits.incl = 10,
                           covFormSite = c(""),
                           covCorGroupSite = c(""),
                           covFormVisit = c("linear"),
                           covCorGroupVisit = c("A"))
  
  # Set up models ----
  
  for (mod in 1:6) {
    
    if (mod == 1) {
      ## PO data ----
      
      # format data
      datalist <- list(PO = po_dat)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simPO1"),
                              data.type = c("PO"),
                              covar.mean = c("POcovA"),
                              covar.sum = c(NA),
                              PO.extent = c("CONUS"))
    }
    
    if (mod == 2) {
      ## count 1 visit data ----
      
      # format data
      datalist <- list(count = count_dat1)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simCount1"),
                              data.type = c("count"),
                              covar.mean = c("cov.siteA"),
                              covar.sum = c(NA),
                              PO.extent = c(NA))
      
      min.visits.incl <- Inf
    }
    
    if (mod == 3) {
      ## count multi visit data ----
      
      # format data
      datalist <- list(count = count_dat2)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simCount1"),
                              data.type = c("count"),
                              covar.mean = c("cov.visitA"),
                              covar.sum = c(NA),
                              PO.extent = c(NA))
      
      min.visits.incl <- 10
    }
    
    if (mod == 4) {
      ## DND 1 visit data ----
      
      # format data
      datalist <- list(DND = dnd_dat1)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simDND1"),
                              data.type = c("DND"),
                              covar.mean = c("cov.siteA"),
                              covar.sum = c(NA),
                              PO.extent = c(NA))
      
      min.visits.incl <- Inf
    }
    
    if (mod == 5) {
      ## count multi visit data ----
      
      # format data
      datalist <- list(DND = dnd_dat2)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simDND1"),
                              data.type = c("DND"),
                              covar.mean = c("cov.visitA"),
                              covar.sum = c(NA),
                              PO.extent = c(NA))
      
      min.visits.incl <- 10
    }
    
    if (mod == 6) {
      ## all data ----
      
      # format data
      datalist <- list(PO = po_dat,
                       count = count_dat1,
                       count = count_dat2,
                       DND = dnd_dat1,
                       DND = dnd_dat2)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simPO1", "simCount2", "simCount3", "simDND4", "simDND5"),
                              data.type = c("PO", "count", "count", "DND", "DND"),
                              covar.mean = c("POcovA", "cov.siteA", "cov.visitA", "cov.siteA", "cov.visitA"),
                              covar.sum = c(NA, NA, NA, NA, NA),
                              PO.extent = c("CONUS", NA, NA, NA, NA))
      
      min.visits.incl <- 10
    }
    
    # Fit models ----
    
    # set up blocks and keys
    spatblocks <- make_CV_blocks(region, rows = 5, cols = 5, k = 3)
    gridkey <- make_gridkey(region, spatblocks, fold.out = "none")
    spatRegion <- NULL
    
    # set up covariates
    covs.lin  <- c("covA", "covB")   # linear terms 
    covs.quad <- c("covC")               # quadratic term
    covs.z    <- c(covs.lin, covs.quad)
    covs.z <- covs.z[!is.na(covs.z) & covs.z != ""]
    
    covs.PO   <- c("POcovA")         # no PO effort covariates
    covs.inat <- NULL                # no iNaturalist data
    
    covar_unscaled <- covar
    
    # scale covariates
    numcols <- which(sapply(covar, is.numeric))
    covar[, numcols] <- sapply(covar[, numcols], scale_this)
    
    # Append squared covariate column
    if (length(covs.quad) > 0 & paste0(covs.quad, collapse = "") != "") {
      for (c in 1:length(covs.quad)) {
        covar[,paste0(covs.quad[c], "2")] <- covar[,covs.quad[c]] * covar[,covs.quad[c]]
        covs.z <- c(covs.z, paste0(covs.quad[c], "2"))
      }
    }
    
    # set up species data
    sp.data <- sppdata_for_nimble(species.data = species.data,
                                  region = region,
                                  file.info = file.info,
                                  covar = covar,
                                  stategrid = NULL,    # no state grid for simulated data
                                  statelines.rm = FALSE,
                                  covs.PO = covs.PO,
                                  covs.inat = covs.inat,
                                  min.visits.incl = min.visits.incl,
                                  occ.mod = TRUE,
                                  nmix.mod = TRUE,
                                  keep.conus.grid.id = gridkey$conus.grid.id[gridkey$group == "train"])
    
    # set up data and constants
    tmp <- data_for_nimble(sp.data = sp.data,
                           covar = covar,
                           covs.z = covs.z,
                           sp.auto = T,      
                           coarse.grid = FALSE,
                           region = region,
                           gridkey = gridkey,
                           spatRegion = NULL)
    data      <- tmp$data
    constants <- tmp$constants
    
    # write nimble code, inits, and params
    code <- nimble_code(data = data,
                        constants = constants,
                        path = tempdir(),
                        sp.auto = T,
                        coarse.grid = FALSE,
                        Bprior = "dnorm(0,1)",
                        block.out = "none",
                        zero_mean = T,
                        tau = 1,
                        min.visits.incl = min.visits.incl,
                        occ.mod = TRUE,
                        nmix.mod = TRUE)
    
    inits <- function(x) {nimble_inits(data = data,
                                       constants = constants,
                                       sp.auto = T,
                                       min.visits.incl = min.visits.incl,
                                       occ.mod = TRUE,
                                       nmix.mod = TRUE,
                                       seed = x)}
    
    params <- nimble_params(data = data,
                            constants = constants,
                            lambda = TRUE,
                            XB = TRUE,
                            sp.auto = T,
                            effort = FALSE)
    
    # fit model
    iter   <- 50000
    thin   <- 5
    burnin <- floor(iter * 0.75)
    
    suppressMessages(source("~/GitHub/species-futures/functions/FXN-nimbleParallel.R"))
    
    samples <- nimbleParallel(code = code,
                              data = data,
                              constants = constants,
                              inits = inits,
                              param = params,
                              iter = iter,
                              burnin = burnin,
                              thin = thin)
    
    # summarize output
    samples <- lapply(samples, get_derived, data = data, project = 0, 
                      coarse.grid = F, spatRegion = NULL)
    
    out <- summarize_samples(samples = samples,
                             data = data,
                             constants = constants,
                             project = 0,
                             coarse.grid = FALSE,
                             block.out = "none",
                             gridkey = gridkey,
                             effort = FALSE,
                             cores = 2L)
    
    # compare output
    names(datalist) <- names(species.data$obs)
    
    par <- sim_compare(out, true = sp, plot = "process")$dat %>%
      mutate(model = mod,
             species = s)
    
    parall <- bind_rows(parall, par)
    # sim_compare(out, true = sp, plot = "process")$plot
  }
  
  
}


parall1 <- parall %>%
  mutate(recovered = case_when(true < hi & true > lo ~ 1,
                               T ~ 0),
         mixed = case_when(rhat < 1.05 ~ 1,
                           T ~ 0))

