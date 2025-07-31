## ---------------------------
## Objective: 
##    - Fit NIMBLE model using setup_BLOCK.rdata file
## 
## Input:
##    - setup_BLOCK.rdata
##    - samples_BLOCK_CHAIN.rds OR samples_BLOCK.rds
##
## Output: 
##    - samples_BLOCK_CHAIN.rds OR samples_BLOCK.rds
##
## ---------------------------

start3 <- Sys.time()


## load packages ----
library(tidyverse, quietly = T)
library(magrittr, quietly = T)
library(sf, quietly = T)
library(SpFut.flexiSDM)



# EDIT THIS SECTION ----
num <- 1
block <- "none"
maxchain <- 3
local <- 1
# ---


# Setup ----
args = commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  num = as.numeric(args[1])
  block = as.numeric(args[2])
  maxchain = as.numeric(args[3])
  local = as.numeric(args[4])
  
  if (local == 0) {
    setwd('/caldera/hovenweep/projects/usgs/ecosystems/eesc/rmummah/iSDM-framework/')
  } 
  
}

# Ensure block is coded correctly
if (block == 4) {
  block <- 'none'
}


mods <- read.csv("code/MVPv1.csv") %>% filter(number %in% num)

sp.code <- mods$sp.code[1]
model <- mods$model[1]


# load output directory, setup, and functions
out.dir = paste0('outputs/',num,'_',sp.code,'_',model,'/')
load(paste0(out.dir,'setup_',block,'.Rdata'))




# Load chains (samples) ----

if (file.exists(paste0(out.dir,'samples_',block,'.rds'))) { # If it was run locally, this file should exist
  samples <- readRDS(paste0(out.dir,'samples_',block,'.rds'))
} else {
  # Load individual sample files and save as a list
  samples <- lapply(1:maxchain, 
                    FUN = function(x) {
                      readRDS(paste0(out.dir,'samples_',block,'_',x,'.rds'))
                    })
}



## Calculate derived quantities ----
# This takes longer with the coarse grid.
samples <- lapply(samples, get_derived, data = data, project = project, 
                  coarse.grid = coarse.grid, spatRegion = spatRegion,
                  pathToProj = 'code/03-species-models/setup-projections.R')



## Summarize chains ----
out <- summarize_samples(samples, 
                         data, 
                         constants, 
                         project = project,
                         coarse.grid = coarse.grid,
                         block.out = block.out,
                         gridkey = gridkey, 
                         spatkey = spatRegion$spatkey, 
                         SLURM = ifelse(local == 1, F, T))

save(out, file = paste0(out.dir, "data", blockname, ".rdata"))



# Plot output ----

# Convergence
ggsave(plot_convergence(out), file = paste0(out.dir, "3_parameters-0_convergence-", block, ".jpg"),
       height = 5, width = 12)

if (block.out == "none") {
  cov.labs <- read.csv("data/covariate-labels.csv")
  
  ### Chains ----
  ggsave(plot_chains(samples, data = data, cov.labs = cov.labs, plot = "B", cutoff = 0), file = paste0(out.dir, "3_parameters-a1_chains-B.jpg"), height = 6, width = 8)
  ggsave(plot_chains(samples, data = data, cov.labs = cov.labs, plot = "alpha", cutoff = 0), file = paste0(out.dir, "3_parameters-b1_chains-alpha.jpg"), height = 6, width = 8)
  
  # Posteriors ----
  ggsave(plot_posteriors(samples, data = data, cov.labs = cov.labs, plot = "B", cutoff = 0), file = paste0(out.dir, "3_parameters-a2_posteriors-B.jpg"), height = 6, width = 8)
  ggsave(plot_posteriors(samples, data = data, cov.labs = cov.labs, plot = "alpha", cutoff = 0), file = paste0(out.dir, "3_parameters-b2_posteriors-alpha.jpg"), height = 6, width = 8)
  
  
  
  ### Parameter estimates ----
  ggsave(plot_pars(out, 
                   plot.type = "full", 
                   plot.group = "process", 
                   title = "Process parameter estimates",
                   cov.labs = cov.labs),
         file = paste0(out.dir, "3_parameters-a3_B.jpg"), height = 6, width = 10)
  
  ggsave(plot_pars(out, 
                   plot.type = "full", 
                   plot.group = "dataset", 
                   title = "Dataset parameter estimates"),
         file = paste0(out.dir, "3_parameters-b3_alpha.jpg"), height = 6, width = 10)
  
  ggsave(plot_pars(out, 
                   plot.type = "full", 
                   plot.group = "observation", 
                   title = "Observation intercept estimates",
                   cov.labs = cov.labs),
         file = paste0(out.dir, "3_parameters-c3_observation.jpg"), height = 6, width = 10)
  
  
  ### Marginal effects ----
  ggsave(plot_effects(data, out, breaks = 0.001), file = paste0(out.dir, "/3_parameters-a4-effects.jpg"), height = 7, width = 10)
  
  
  ### Maps ----
  # Current
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"), 
                         region, 
                         plot = "lambda",
                         out = out)
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-a_lambda.jpg"), height = 7, width = 9)
  
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                         region, 
                         plot = "lambda",
                         out = out, 
                         plot.uncertainty = "unc.range")
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-a_lambda-uncabs.jpg"), height = 7, width = 9)
  
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                         region,
                         plot = "lambda",
                         out = out, 
                         plot.uncertainty = "unc.rel")
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-a_lambda-uncrel.jpg"), height = 7, width = 9)
  
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                         region,
                         plot = "psi",
                         out = out)
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-b_psi.jpg"), height = 7, width = 9)
  
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                         region,
                         plot = "psi",
                         out = out, 
                         plot.uncertainty = "unc.range")
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-b_psi-uncabs.jpg"), height = 7, width = 9)
  
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                         region,
                         plot = "psi",
                         out = out, 
                         plot.uncertainty = "unc.rel")
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-b_psi-uncrel.jpg"), height = 7, width = 9)
  
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                         region,
                         plot = "boundary",
                         out = out, 
                         threshold = 0.25)
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-c1_boundary-0.25.jpg"), height = 7, width = 9)
  
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                         region,
                         plot = "boundary",
                         out = out, 
                         threshold = 0.5)
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-c2_boundary-0.5.jpg"), height = 7, width = 9)
  
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                         region,
                         plot = "boundary",
                         out = out, 
                         threshold = 0.75)
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-c3_boundary-0.75.jpg"), height = 7, width = 9)
  
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                         region,
                         plot = "spat",
                         out = out,
                         coarse.grid = coarse.grid,
                         spatRegion = spatRegion)
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-d_spat.jpg"), height = 7, width = 9)
  
  pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                         region,
                         plot = "XB",
                         out = out,
                         plot.exp = T)
  ggsave(pl, file = paste0(out.dir, "4_mapcurrent-d_expXB.jpg"), height = 7, width = 9)
  
  
  
  ### Future projections ----
  if (project > 0) {
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "lambda",
                           out = out, 
                           plot.current = F, 
                           plot.change = F,
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "5_mapfuture-a_lambda.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "lambda",
                           out = out, 
                           plot.current = F, 
                           plot.change = F, 
                           plot.uncertainty = "unc.range",
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "5_mapfuture-a_lambda-uncabs.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "lambda",
                           out = out, 
                           plot.current = F, 
                           plot.change = F, 
                           plot.uncertainty = "unc.rel",
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "5_mapfuture-a_lambda-uncrel.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region, 
                           plot = "psi",
                           out = out, 
                           plot.current = F, 
                           plot.change = F,
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "5_mapfuture-b_psi.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "psi",
                           out = out, 
                           plot.current = F, 
                           plot.change = F, 
                           plot.uncertainty = "unc.range",
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "5_mapfuture-b_psi-uncabs.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "psi",
                           out = out, 
                           plot.current = F, 
                           plot.change = F, 
                           plot.uncertainty = "unc.rel",
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "5_mapfuture-b_psi-uncrel.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "boundary",
                           out = out, 
                           plot.current = F, 
                           plot.change = F, 
                           threshold = 0.25,
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "5_mapfuture-c1_boundary-0.25.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "boundary",
                           out = out, 
                           plot.current = F, 
                           plot.change = F, 
                           threshold = 0.5,
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "5_mapfuture-c2_boundary-0.5.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "boundary",
                           out = out, 
                           plot.current = F, 
                           plot.change = F, 
                           threshold = 0.75,
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "5_mapfuture-c3_boundary-0.75.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "XB",
                           out = out, 
                           plot.exp = T,
                           plot.current = F,
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "5_mapfuture-d_expXB-future.jpg"), height = 7, width = 9)
    
    
    ### Future absolute change ----
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "lambda",
                           out = out, 
                           plot.current = F, 
                           plot.change = "absolute",
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "6_mapchange-a1_lambda-abs.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"),
                           region,
                           plot = "psi",
                           out = out, 
                           plot.current = F, 
                           plot.change = "absolute",
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "6_mapchange-b1_psi-abs.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"), 
                           region,
                           plot = "boundary",
                           out = out, 
                           plot.current = F, 
                           plot.change = "absolute", 
                           threshold = 0.25,
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "6_mapchange-c1_boundary-abs-0.25.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"), 
                           region,
                           plot = "boundary",
                           out = out, 
                           plot.current = F, 
                           plot.change = "absolute", 
                           threshold = 0.5,
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "6_mapchange-c2_boundary-abs-0.5.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"), 
                           region,
                           plot = "boundary",
                           out = out, 
                           plot.current = F, 
                           plot.change = "absolute", 
                           threshold = 0.75,
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "6_mapchange-c3_boundary-abs-0.75.jpg"), height = 7, width = 9)
    
    ### Future relative change ----
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"), 
                           region,
                           plot = "lambda",
                           out = out, 
                           plot.current = F, 
                           plot.change = "relative",
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "6_mapchange-a2_lambda-rel.jpg"), height = 7, width = 9)
    
    pl <- map_species_data(title = paste0(common, " (", sp.code, ")"), 
                           region,
                           plot = "psi", 
                           out = out, 
                           plot.current = F,
                           plot.change = "relative",
                           proj.names = proj.names)
    ggsave(pl, file = paste0(out.dir, "6_mapchange-b2_psi-rel.jpg"), height = 7, width = 9)
  }
}



# Validate ----

# Get all data
species.data <- load_species_data(sp.code = sp.code,
                                  sp.code.all = sp.code.all,
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



val.dat <- c()
for (d in 1:length(species.data$obs)) {
  tmp <- species.data$obs[[d]] %>%
    mutate(site.id = as.character(site.id))
  val.dat <- bind_rows(val.dat, tmp)
}

val.dat <- val.dat %>%
  filter(data.type != "PO") %>%
  mutate(pa = case_when(count > 0 ~ 1,
                        count == 0 ~ 0)) %>%
  select(conus.grid.id, source, pa)


# If it's only PO data, skip the following
if (nrow(val.dat) == 0) {
  
  # Calculate AUC for in sample
  val.in <- filter(out$psi0, group == "train") %>%
    inner_join(val.dat, by = "conus.grid.id")
  
  # overall
  if (length(unique(val.in$pa)) < 2) {
    AUCin <- NA
  } else {
    AUCin <- pROC::auc(val.in$pa, val.in$mean)
  }
  
  # each dataset
  dats <- unique(val.in$source)
  aucsin <- c()
  for (s in 1:length(dats)) {
    val.in1 <- filter(val.in, source == dats[s])
    
    if (length(unique(val.in1$pa)) < 2) {
      AUCin1 <- NA
    } else {
      AUCin1 <- pROC::auc(val.in1$pa, val.in1$mean)
    }
    
    tmpin <- data.frame(source = dats[s],
                        AUCin = as.numeric(AUCin1),
                        in.n = nrow(val.in1),
                        in.cell = length(unique(val.in1$conus.grid.id)))
    
    aucsin <- bind_rows(aucsin, tmpin)
  }
  
  
  
  # Calculate AUC for out of sample
  val.out <- filter(out$psi0, group == "test") %>% inner_join(val.dat, by = "conus.grid.id")
  
  # Overall
  if (nrow(val.out) == 0) {
    # if there is no out of sample data
    AUCout <- NA
    
  } else {
    # if there is out of sample data
    AUCout <- pROC::auc(val.out$pa, val.out$mean)
  }

  
  # each dataset
  dats <- unique(val.out$source)
  aucsout <- c()
  for (s in 1:length(dats)) {
    val.out1 <- filter(val.out, source == dats[s])
    
    if (length(unique(val.out1$pa)) < 2) {
      AUCout1 <- NA
    } else {
      AUCout1 <- pROC::auc(val.out1$pa, val.out1$mean)
    }
    
    tmpout <- data.frame(source = dats[s],
                         AUCout = as.numeric(AUCout1),
                         out.n = nrow(val.out1),
                         out.cell = length(unique(val.out1$conus.grid.id)))
    
    aucsout <- bind_rows(aucsout, tmpout)
    
  }
  
  
  
  if (block.out != "none") {aucs <- full_join(aucsin, aucsout, by = "source")}
  if (block.out == "none") {aucs <- aucsin}
  
  all.auc <- data.frame(sp.code = sp.code,
                        block = block.out,
                        AUCin.full = as.numeric(AUCin), 
                        in.full.n = nrow(val.in),
                        in.full.cell = length(unique(val.in$conus.grid.id)),
                        AUCout.full = as.numeric(AUCout),
                        out.full.n = nrow(val.out),
                        out.full.cell = length(unique(val.out$conus.grid.id)),
                        aucs)
  
}

# If there are no validation data, set auc.out = NA
if (!exists("all.auc")) {
  all.auc <- NA
}


# Save everything ----
if (block.out == "none") {
  blockname <- "full"
  save(region, file = paste0(out.dir, "region.rdata"))
} else {blockname <- block.out}


save(species.data, covar, covar_unscaled, data, constants, gridkey, all.auc, block.out,
     file = paste0(out.dir, "data", blockname, "-info.rdata"))



# Final CV Figures ----
# Make final cross validation figs if all blocks have been done

done <- list.files(path = out.dir)

mod.name <- c('data1','data2','data3','datafull')


if (all(c(paste0(mod.name, '.rdata'), paste0(mod.name, '-info.rdata')) %in% done)) {
  
  
  # Load data from all model runs
  auc <- c()
  proc <- c()
  obs <- c()
  alpha <- c()
  for (i in 1:3) {
    load(paste0(out.dir, "data", i, "-info.rdata"))
    load(paste0(out.dir, "data", i, ".rdata"))
    
    # If validation data is missing, skip
    if (is.na(all.auc)) next
    
    out$process.coef$block.out <- as.character(out$process.coef$block.out)
    out$obs.coef$block.out <- as.character(out$obs.coef$block.out)
    out$alpha$block.out <- as.character(out$alpha$block.out)
    
    auc <- bind_rows(auc, all.auc)
    proc <- bind_rows(proc, out$process.coef)
    obs <- bind_rows(obs, out$obs.coef)
    alpha <- bind_rows(alpha, out$alpha)
    
  }
  
  # now add full model
  load(paste0(out.dir, "datafull.rdata"))
  
  out$process.coef$block.out <- as.character(out$process.coef$block.out)
  out$obs.coef$block.out <- as.character(out$obs.coef$block.out)
  out$alpha$block.out <- as.character(out$alpha$block.out)
  
  proc <- bind_rows(proc, out$process.coef)
  obs <- bind_rows(obs, out$obs.coef)
  alpha <- bind_rows(alpha, out$alpha)
  
  
  load(paste0(out.dir, "datafull-info.rdata"))
  auc$block <- as.character(auc$block)
  auc <- bind_rows(auc, all.auc)
  
  
  ggsave(plot_pars(out = proc, 
                   plot.type = "cv", 
                   plot.group = "process", 
                   title = "Process parameter estimates",
                   cov.labs = cov.labs),
         file = paste0(out.dir, "7_CV-a-parameters-process.jpg"), height = 6, width = 10)
  
  ggsave(plot_pars(out = alpha, 
                   plot.type = "cv", 
                   plot.group = "dataset", 
                   title = "Dataset parameter estimates",
                   cov.labs = cov.labs),
         file = paste0(out.dir, "7_CV-b-parameters-dataset.jpg"), height = 6, width = 10)
  
  ggsave(plot_pars(out = obs, 
                   plot.type = "cv", 
                   plot.group = "observation", 
                   title = "Observation parameter estimates",
                   cov.labs = cov.labs),
         file = paste0(out.dir, "7_CV-a-parameters-observation.jpg"), height = 6, width = 10)
  
  
  # Plot AUC
  ggsave(plot_auc(auc, lines = T), file = paste0(out.dir, "7_CV-d_AUClines.jpg"), height = 7, width = 10)
  ggsave(plot_auc(auc, lines = F), file = paste0(out.dir, "7_CV-d_AUCnolines.jpg"), height = 7, width = 10)
  
}


# Create time-summary output file ----
cat("Summarization and visualization complete \n\n")

end3 <- Sys.time() - start3

# Start output file
sink(paste0(out.dir,'time-summary_',block,'.txt'))

cat(paste0("Time summary for: ", nums.do,' ',sp.code,' ',model,'\n\n'))

cat(paste0('The model was run on ', format(Sys.time(), "%a %b %Y %d %X"),'\n\n'))


# Print total time to output file
cat('Total time for setup:\n')
print(end1)
cat('\n')

end2 <- readRDS(paste0(out.dir, 'time2_',block,'.rds'))
cat('Total time for fitting NIMBLE model: \n')
print(end2)
cat('\n')

cat('Total time for summarization: \n')
print(end3)
cat('\n')

cat('Start to finish time (code): \n')
print(end1 + end2 + end3)
cat('\n')

cat('Start to finish time (user delay): \n')
print(Sys.time() - start1)
cat('\n\n')

# Close output file
sink()

cat("Removing intermediate files...\n")
# system(paste0("rm ./outputs/03-species-models/MVPv1/",number,'_',sp.code,'_',model,"/setup_",block,"*"))
# system(paste0("rm ./outputs/03-species-models/MVPv1/",number,'_',sp.code,'_',model,"/samples_",block,"*"))
# system(paste0("rm ./outputs/03-species-models/MVPv1/",number,'_',sp.code,'_',model,"/time2_",block,".rds"))

cat(paste0("Model ",number,' ',sp.code,' ',model," complete!!\n"))


# End script

