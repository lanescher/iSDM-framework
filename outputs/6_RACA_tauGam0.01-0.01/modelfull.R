code <- nimbleCode({


# Process Model

# 21560 cells
for (i in 1:nCell) {

    # log intensity
    XB0[i] <- inprod(B[1:nCovZ], Xz[i,1:nCovZ])
    log(lambda0[i]) <- XB0[i] + spat[i] # residual spatial effect

}


# ---------------------------------------

# Observation Model 1: PO, BLM Surveys
# 21560 cells
for (j in 1:nW1) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD1[j] <- lambda0[Wcells1[j]] * alpha[1]

  # Observation model
  W1[j] ~ dpois(lambdaD1[j] * E1[j]) # Poisson

  # Effort 
  log(E1[j]) <- A1[1] * Xw1[j,1]

  # Prior for X imputation
  Xw1[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 2: PO, NCCN
# 21560 cells
for (j in 1:nW2) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD2[j] <- lambda0[Wcells2[j]] * alpha[2]

  # Observation model
  W2[j] ~ dpois(lambdaD2[j] * E2[j]) # Poisson

  # Effort 
  log(E2[j]) <- A2[1] * Xw2[j,1]

  # Prior for X imputation
  Xw2[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 3: PO, USGS RACA assessment - current, FS NRIS, ORBIC, WDWF
# 21560 cells
for (j in 1:nW3) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD3[j] <- lambda0[Wcells3[j]] * alpha[3]

  # Observation model
  W3[j] ~ dpois(lambdaD3[j] * E3[j]) # Poisson

  # Effort 
  log(eff3[j]) <- inprod(A3[1:nCovW3], Xw3[j,1:nCovW3])
  E3[j] <- eff3[j] * S3[j]
}

  

# Observation Model 4: DND, ARMI VES
# 800 observations, 5 median visits per site
for (j in 1:nV4) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD4[j] <- lambda0[Vcells4[j]] * alpha[4]

  # Observation model
  V4[j] ~ dbern(ZD4[j] * p4[j]) # Occupancy

  # Detection 
  logit(p4[j]) <- inprod(D4[1:nCovV4], Xv4[j,1:nCovV4])

  # Prior for X imputation
  for (c in 1:nCovV4) {
    Xv4[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 5: count, PNW Salamanders
# 76 observations, 2 median visits per site
for (j in 1:nY5) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD5[j] <- lambda0[Ycells5[j]] * alpha[5]

  # Observation model
  Y5[j] ~ dpois(lambdaD5[j] * p5) # Poisson

  
}
# Detection 
log(p5) <- 0
  

# Observation Model 6: count, USGS RACA assessment - historic
# 115 observations, 2 median visits per site
for (j in 1:nY6) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD6[j] <- lambda0[Ycells6[j]] * alpha[6]

  # Observation model
  Y6[j] ~ dpois(lambdaD6[j] * p6) # Poisson

  
}
# Detection 
log(p6) <- 0
  
# ---------------------------------------
# Process priors
for (a in 1:nCovZ) {
  B[a] ~ dnorm(0,1)
}

# Dataset intercepts
for (a in 1:nD) {
  alpha[a] <- exp(w[a])
  w[a] ~ dnorm(0,1)
}

# Observation priors, PO 1: BLM Surveys
for (b in 1:nCovW1) {
  A1[b] ~ dnorm(0,1)
}

# Observation priors, PO 2: NCCN
for (b in 1:nCovW2) {
  A2[b] ~ dnorm(0,1)
}

# Observation priors, PO 3: USGS RACA assessment - current, FS NRIS, ORBIC, WDWF
for (b in 1:nCovW3) {
  A3[b] ~ dnorm(0,1)
}

# Observation priors, DND 4: ARMI VES
for (b in 1:nCovV4) {
  D4[b] ~ dnorm(0,1)
}



tau ~ dgamma(0.01,0.01)
spat[1:nCell] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nCell], tau, zero_mean = 1)


})
