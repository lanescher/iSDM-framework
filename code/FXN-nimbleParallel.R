
# https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html

run_nimbleMCMC <- function(info, code, constants, data, param, 
                           iter, burnin, thin) {
  
  library(nimble)
  # library(nimbleHMC)
  
  # samples <- nimbleMCMC(code = code,
  #                       constants = constants,
  #                       data = data,
  #                       inits = info$inits, 
  #                       monitors = param,
  #                       niter = iter,
  #                       nburnin = burnin,
  #                       nchains = 1,
  #                       thin = thin,
  #                       setSeed = info$seed,
  #                       samplesAsCodaMCMC = TRUE)
  
  Rmodel <- nimbleModel(code = code, 
                        name = 'Rmodel',
                        constants = constants, 
                        data = data,
                        inits = info$inits, 
                        check = FALSE) 
  # Must have this for HMC
  # buildDerivs = T)
  # Rmodel$initializeInfo()
  conf <- configureMCMC(Rmodel, monitors = param)
  # conf$addSampler(target = 'B', type = 'NUTS')
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  ## Alternative way of compiling/building after the Rmodel step
  # Cmodel <- compileNimble(Rmodel)
  
  # mcmc <- buildMCMC(Cmodel, monitors = param)
  
  # Cmcmc <- compileNimble(mcmc)
  
  samples <- runMCMC(Cmcmc, 
                     niter = iter,  
                     nburnin = burnin, 
                     thin = thin,
                     setSeed = info$seed,
                     samplesAsCodaMCMC = TRUE)
  
  return(samples)
}


# removed core argument
nimbleParallel <- function(cores = 3, code, constants, data, 
                           inits, param, iter, burnin, thin) {
  library(parallel)

  this_cluster <- makeCluster(cores)
  
  tmp <- lapply(1:cores, function(x) list(seed = x,
                                          inits = inits(x)))
  
  # run in parallel
  cat("Jack be NIMBLE... \n")
  cat("Jack be quick... \n")
  
  chain_output <- parLapply(cl = this_cluster,
                            X = tmp,
                            fun = run_nimbleMCMC,
                            code = code,
                            constants = constants,
                            data = data,
                            param = param,
                            iter = iter,
                            burnin = burnin,
                            thin = thin)
  
  cat("Jack jumped over the candlestick! \n")
  
  stopCluster(this_cluster)
  
  return(chain_output)
}