## ---------------------------
## Objective: 
##    - Fit NIMBLE model using setup_BLOCK.rdata file
## 
## Input:
##    - setup_BLOCK.rdata
##
## Output: 
##    - samples_BLOCK_CHAIN.rds OR samples_BLOCK.rds
##
## ---------------------------

start2 <- Sys.time() 
print(paste0('Beginning 02-flexiSDM script at ', start2))


# EDIT THIS SECTION ----
num <- 3
block <- "none"
chain <- 1
local <- 1
# ---


args = commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  num = as.numeric(args[1])
  block = as.numeric(args[2])
  chain = as.numeric(args[3])
  local = as.numeric(args[4]) 
  
  # If running on HPC, set working directory
  if (local == 0) {
    setwd('/caldera/hovenweep/projects/usgs/ecosystems/eesc/cscher/iSDM-framework/')
  } 
}

# Ensure block is coded correctly
if (block == 4) {
  block <- 'none'
}


# Load functions and packages
suppressMessages(source("code/FXN-nimbleParallel.R"))
library(SpFut.flexiSDM)
library(tidyverse)


mods <- read.csv("code/model-specs.csv") %>% filter(number %in% num)

sp.code <- mods$sp.code[1]
model <- mods$model[1]

# Set output directory and load setup file
out.dir = paste0('outputs/',num,'_',sp.code,'_',model,'/')
load(paste0(out.dir,'setup_',block,'.Rdata'))


# Fit NIMBLE model ----
if (local == 1) {
  start.nim <- Sys.time()
  samples <- nimbleParallel(code = code,
                            data = data,
                            constants = constants,
                            inits = inits,
                            param = params,
                            iter = 5000,
                            burnin = 1000,
                            thin = thin)
  
  end.nim <- Sys.time() - start.nim
  print(end.nim)
  
  saveRDS(samples, paste0(out.dir,'samples_',block,'.rds'))
} else {

  info <- list(seed = chain,
               inits = inits(chain))
  
  start.nim <- Sys.time()
  samples <- run_nimbleMCMC(info = info, 
                            code = code, 
                            constants = constants, 
                            data = data, 
                            param = params, 
                            iter = iter, 
                            burnin = burnin, 
                            thin = thin)
  
  end.nim <- Sys.time() - start.nim
  print(end.nim)
  
  saveRDS(samples, paste0(out.dir,'samples_',block,'_',chain,'.rds'))
}


end2 <- Sys.time() - start2
saveRDS(end2, paste0(out.dir,'time2_',block,'.rds'))
print(end2)

# End script - proceed to 03-flexiSDM.R
