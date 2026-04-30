

# Simulation study for model framework paper


# remotes::install_github("rileymummah/SpFut.flexiSDM-dev@v1.1-updates", force = T, auth_token = Sys.getenv("pat"))
library(SpFut.flexiSDM)
library(sf)
library(nimble)
library(tidyverse)



# parameters for species
grid.dim <- 20
shape <- "hexagon"
covForm <- c("linear", "linear", "quadratic")
covCorGroup <- c("A", "B", "C")
spatial <- T



parall <- c()
for (s in 1) {
  
  set.seed(s)
  
  # Simulate species ----
  sp <- sim_species(grid.dim = grid.dim, 
                    shape = shape,
                    cellsize.sim = 1,         
                    cellsize.inf = 1,        
                    covForm = covForm,
                    covCorGroup = covCorGroup,
                    covType = rep("spatial", times = length(covForm)),
                    covNoise = rep(NA, times = length(covForm)),
                    paramShift = 0,
                    paramSD = 0.1,
                    spatial = spatial)
  
  # Simulate data ----
  
  ## PO data ----
  po_dat <- sim_POdat(true = sp, 
                      shape = shape,
                      link = "log",
                      covForm = c("linear"),
                      covCorGroup = c("A"),
                      covType = c("spatial"),
                      covNoise = c(NA),
                      paramSD = 0.5)
  
  ## count 1 visit - no covs ----
  count_dat1a <- sim_surveydat(true = sp,
                             survey.type = "count",
                             method = "random",
                             nsite = 1000,
                             nvisit = 1,
                             nvisit.sd = 0,
                             link = "log",
                             min.visits.incl = Inf,
                             covFormSite = c(""),
                             covCorGroupSite = c(""),
                             covFormVisit = c(""),
                             covCorGroupVisit = c(""))
  
  ## count 1 visit ----
  count_dat1b <- sim_surveydat(true = sp,
                              survey.type = "count",
                              method = "random",
                              nsite = 1000,
                              nvisit = 1,
                              nvisit.sd = 0,
                              link = "log",
                              min.visits.incl = Inf,
                              covFormSite = c("linear"),
                              covCorGroupSite = c("A"),
                              covFormVisit = c(""),
                              covCorGroupVisit = c(""))
  
  ## count 10 visits - no covs ----
  count_dat2a <- sim_surveydat(true = sp,
                              survey.type = "count",
                              method = "random",
                              nsite = 1000,
                              nvisit = 10,
                              nvisit.sd = 0,
                              link = "logit",
                              min.visits.incl = 3,
                              covFormSite = c(""),
                              covCorGroupSite = c(""),
                              covFormVisit = c(""),
                              covCorGroupVisit = c(""))
  
  ## count 10 visits ----
  count_dat2b <- sim_surveydat(true = sp,
                             survey.type = "count",
                             method = "random",
                             nsite = 1000,
                             nvisit = 10,
                             nvisit.sd = 0,
                             link = "logit",
                             min.visits.incl = 3,
                             covFormSite = c("linear"),
                             covCorGroupSite = c("A"),
                             covFormVisit = c(""),
                             covCorGroupVisit = c(""))
  
  ## DND 1 visit - no covs ----
  dnd_dat1a <- sim_surveydat(true = sp,
                           survey.type = "DND",
                           method = "regular",
                           nsite = 1000,
                           nvisit = 1,
                           nvisit.sd = 0,
                           link = "log",
                           min.visits.incl = Inf,
                           covFormSite = c(""),
                           covCorGroupSite = c(""),
                           covFormVisit = c(""),
                           covCorGroupVisit = c(""))
  
  ## DND 1 visit ----
  dnd_dat1b <- sim_surveydat(true = sp,
                            survey.type = "DND",
                            method = "regular",
                            nsite = 1000,
                            nvisit = 1,
                            nvisit.sd = 0,
                            link = "log",
                            min.visits.incl = Inf,
                            covFormSite = c("linear"),
                            covCorGroupSite = c("A"),
                            covFormVisit = c(""),
                            covCorGroupVisit = c(""))
  
  
  ## DND 10 visits - no covs ----
  dnd_dat2a <- sim_surveydat(true = sp,
                           survey.type = "DND",
                           method = "regular",
                           nsite = 1000,
                           nvisit = 10,
                           nvisit.sd = 0,
                           link = "logit",
                           min.visits.incl = 3,
                           covFormSite = c(""),
                           covCorGroupSite = c(""),
                           covFormVisit = c(""),
                           covCorGroupVisit = c(""))
  
  
  ## DND 10 visits ----
  dnd_dat2b <- sim_surveydat(true = sp,
                            survey.type = "DND",
                            method = "regular",
                            nsite = 1000,
                            nvisit = 10,
                            nvisit.sd = 0,
                            link = "logit",
                            min.visits.incl = 3,
                            covFormSite = c("linear"),
                            covCorGroupSite = c("A"),
                            covFormVisit = c(""),
                            covCorGroupVisit = c(""))
  
  # Set up models ----
  
  for (mod in 1:11) {
    
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
      ## count 1 visit data, no covs ----
      
      # format data
      datalist <- list(count = count_dat1a)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simCount1"),
                              data.type = c("count"),
                              covar.mean = c(""),
                              covar.sum = c(NA),
                              PO.extent = c(NA))
      
      min.visits.incl <- Inf
    }
    
    if (mod == 3) {
      ## count 1 visit data ----
      
      # format data
      datalist <- list(count = count_dat1b)
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
    
    if (mod == 4) {
      ## count multi visit data, no covs ----
      
      # format data
      datalist <- list(count = count_dat2a)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simCount1"),
                              data.type = c("count"),
                              covar.mean = c(""),
                              covar.sum = c(NA),
                              PO.extent = c(NA))
      
      min.visits.incl <- 3
    }
    
    if (mod == 5) {
      ## count multi visit data ----
      
      # format data
      datalist <- list(count = count_dat2b)
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
      
      min.visits.incl <- 3
    }
    
    if (mod == 6) {
      ## DND 1 visit data, no covs ----
      
      # format data
      datalist <- list(DND = dnd_dat1a)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simDND1"),
                              data.type = c("DND"),
                              covar.mean = c(""),
                              covar.sum = c(NA),
                              PO.extent = c(NA))
      
      min.visits.incl <- Inf
    }
    
    if (mod == 7) {
      ## DND 1 visit data ----
      
      # format data
      datalist <- list(DND = dnd_dat1b)
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
    
    if (mod == 8) {
      ## count multi visit data, no covs ----
      
      # format data
      datalist <- list(DND = dnd_dat2a)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simDND1"),
                              data.type = c("DND"),
                              covar.mean = c(""),
                              covar.sum = c(NA),
                              PO.extent = c(NA))
      
      min.visits.incl <- 3
    }
    
    if (mod == 9) {
      ## count multi visit data ----
      
      # format data
      datalist <- list(DND = dnd_dat2b)
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
      
      min.visits.incl <- 3
    }
    
    if (mod == 10) {
      ## all data, no covs ----
      
      # format data
      datalist <- list(PO = po_dat,
                       count = count_dat1a,
                       count = count_dat2a,
                       DND = dnd_dat1a,
                       DND = dnd_dat2a)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simPO1", "simCount2", "simCount3", "simDND4", "simDND5"),
                              data.type = c("PO", "count", "count", "DND", "DND"),
                              covar.mean = c("POcovA", "", "", "", ""),
                              covar.sum = c(NA, NA, NA, NA, NA),
                              PO.extent = c("CONUS", NA, NA, NA, NA))
      
      min.visits.incl <- 3
    }
    
    if (mod == 11) {
      ## all data ----
      
      # format data
      datalist <- list(PO = po_dat,
                       count = count_dat1b,
                       count = count_dat2b,
                       DND = dnd_dat1b,
                       DND = dnd_dat2b)
      sim_out <- format_sim(true = sp,
                            data = datalist)
      
      region       <- sim_out$region
      species.data <- sim_out$species.data
      covar        <- sim_out$covar
      
      file.info <- data.frame(file.label = c("simPO1", "simCount2", "simCount3", "simDND4", "simDND5"),
                              data.type = c("PO", "count", "count", "DND", "DND"),
                              covar.mean = c("POcovA", "cov.siteA", "cov.siteA", "cov.siteA", "cov.siteA"),
                              covar.sum = c(NA, NA, NA, NA, NA),
                              PO.extent = c("CONUS", NA, NA, NA, NA))
      
      min.visits.incl <- 3
    }
    
    # Fit models ----
    
    # set up blocks and keys
    spatblocks <- make_CV_blocks(region, rows = 5, cols = 5, k = 3)
    gridkey <- make_gridkey(region, spatblocks, fold.out = "none")
    spatRegion <- NULL
    
    # set up covariates
    covs.linear <- which(covForm == "linear")
    covs.quadratic <- which(covForm == "quadratic")
    covs.lin  <- colnames(sp$covar)[covs.linear + 1]   # linear terms 
    covs.quad <- colnames(sp$covar)[covs.quadratic + 1]
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
                           sp.auto = spatial,      
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
                        sp.auto = spatial,
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
                                       sp.auto = spatial,
                                       min.visits.incl = min.visits.incl,
                                       occ.mod = TRUE,
                                       nmix.mod = TRUE,
                                       seed = x)}
    
    params <- nimble_params(data = data,
                            constants = constants,
                            lambda = TRUE,
                            XB = TRUE,
                            sp.auto = spatial,
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
             species = s,
             type = "process")
    parall <- bind_rows(parall, par)
    
    par <- sim_compare(out, true = sp, plot = "alpha")$dat %>%
      mutate(model = mod,
             species = s,
             type = "alpha")
    parall <- bind_rows(parall, par)
    
    # par <- sim_compare(out, true = sp, plot = "obs")$dat %>%
    #   mutate(model = mod,
    #          species = s,
    #          type = "obs")
    # parall <- bind_rows(parall, par)
    
    
    # sim_compare(out, true = sp, plot = "process")$plot

  }
  
  
}


parall1 <- parall %>%
  mutate(recovered = case_when(true < hi & true > lo ~ 1,
                               T ~ 0),
         mixed = case_when(rhat < 1.05 ~ 1,
                           T ~ 0),
         model.lab = case_when(model == 1 ~ "PO",
                               model == 2 ~ "count, 1 visit",
                               model == 3 ~ "count, 1 visit, covs",
                               model == 4 ~ "count, 10 visits",
                               model == 5 ~ "count, 10 visits, covs",
                               model == 6 ~ "DND, 1 visit",
                               model == 7 ~ "DND, 1 visit, covs",
                               model == 8 ~ "DND, 10 visits",
                               model == 9 ~ "DND, 10 visits, covs",
                               model == 10 ~ "all", 
                               model == 11 ~ "all, covs"))


ggplot(filter(parall1, type == "process", mixed == 1)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x = covariate, y = true), size = 3, shape = 9, stroke = 1) +
  geom_pointrange(aes(x = covariate, y = mean, ymin = lo, ymax = hi, color = model.lab), position = position_dodge(width = 0.5)) +
  facet_wrap(~ species) +
  theme_bw() +
  labs(x = "Covariate", y = "Estimate", color = "Model") +
  coord_cartesian(ylim = c(-0.5, 0.5))


cov.labs <- data.frame(covariate = c(covs.z, covs.PO, covs.inat),
                       Label = c(covs.z, covs.PO, covs.inat))

outplot <- plot_chains(samples, data = data, cov.labs = cov.labs,
                       plot = "B", cutoff = 0)
outplot$plot
