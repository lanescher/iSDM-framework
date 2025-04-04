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
print(paste0('Beginning 02-MVPv1.1 script at ', start2))


# EDIT THIS SECTION ----
num <- 1
block <- "none"
sp.code <- "RACA"
model <- "WIPtest"
chain <- 1
local <- 1
# ---


args = commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  num = as.numeric(args[1])
  sp.code = as.character(args[2])
  model = as.character(args[3])
  block = as.numeric(args[4])
  chain = as.numeric(args[5])
  local = as.numeric(args[6]) 
  
  # If running on HPC, set working directory
  if (local == 0) {
    setwd('/caldera/hovenweep/projects/usgs/ecosystems/eesc/rmummah/species-futures/')
  } 
}

# Ensure block is coded correctly
if (block == 4) {
  block <- 'none'
}


# Load functions and packages
suppressMessages(source("functions/FXN-nimbleParallel.R"))
library(SpFut.flexiSDM)


# Set output directory and load setup file
out.dir = paste0('outputs/03-species-models/MVPv1/',num,'_',sp.code,'_',model,'/')
load(paste0(out.dir,'setup_',block,'.Rdata'))


# Fit NIMBLE model ----
if (local == 1) {
  start.nim <- Sys.time()
  samples <- nimbleParallel(code = code,
                            data = data,
                            constants = constants,
                            inits = inits,
                            param = params,
                            iter = iter,
                            burnin = burnin,
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