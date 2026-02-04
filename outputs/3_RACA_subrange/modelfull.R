code <- nimble::nimbleCode({


# Process Model

# 279 cells
for (i in 1:nCell) {

    # log intensity
    XB0[i] <- inprod(B[1:nCovZ], Xz[i,1:nCovZ])
    log(lambda0[i]) <- XB0[i] + spat[i] # residual spatial effect

}


# ---------------------------------------

# Observation Model 1: PO, FS NRIS
# 279 cells
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

  

# Observation Model 2: PO, BLM Surveys, ORBIC
# 279 cells
for (j in 1:nW2) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD2[j] <- lambda0[Wcells2[j]] * alpha[2]

  # Observation model
  W2[j] ~ dpois(lambdaD2[j] * E2[j]) # Poisson

  # Effort 
  log(eff2[j]) <- A2[1] * Xw2[j,1]
  E2[j] <- eff2[j] * S2[j]

  # Prior for X imputation
  Xw2[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 3: DND, USGS RACA assessment - historic
# 24 observations, 2.5 median visits per site
for (j in 1:nV3) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD3[j] <- lambda0[Vcells3[j]] * alpha[3]

  # Observation model
  V3[j] ~ dbern(1-exp(-lambdaD3[j] * p3[j])) # Bernoulli

  # Detection 
  log(p3[j]) <- inprod(D3[1:nCovV3], Xv3[j,1:nCovV3])

  # Prior for X imputation
  for (c in 1:nCovV3) {
    Xv3[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 4: count, USGS RACA assessment - current
# 16 observations, 2 median visits per site
for (j in 1:nY4) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD4[j] <- lambda0[Ycells4[j]] * alpha[4]

  # Observation model
  Y4[j] ~ dpois(lambdaD4[j] * p4[j]) # Poisson

  # Detection 
  log(p4[j]) <- C4[1] * Xy4[j,1]

  # Prior for X imputation
  Xy4[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 5: count, PNW Salamanders
# 19 observations, 1 median visits per site
for (j in 1:nY5) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD5[j] <- lambda0[Ycells5[j]] * alpha[5]

  # Observation model
  Y5[j] ~ dpois(lambdaD5[j] * p5[j]) # Poisson

  # Detection 
  log(p5[j]) <- C5[1] * Xy5[j,1]

  # Prior for X imputation
  Xy5[j, 1] ~ dnorm(0, 1)
}

  
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

# Observation priors, PO 1: FS NRIS
for (b in 1:nCovW1) {
  A1[b] ~ dnorm(0,1)
}

# Observation priors, PO 2: BLM Surveys, ORBIC
for (b in 1:nCovW2) {
  A2[b] ~ dnorm(0,1)
}

# Observation priors, DND 3: USGS RACA assessment - historic
for (b in 1:nCovV3) {
  D3[b] ~ dnorm(0,1)
}

# Observation priors, count 4: USGS RACA assessment - current
for (b in 1:nCovY4) {
  C4[b] ~ dnorm(0,1)
}

# Observation priors, count 5: PNW Salamanders
for (b in 1:nCovY5) {
  C5[b] ~ dnorm(0,1)
}

tau <- 1
spat[1:nCell] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nCell], tau, zero_mean = 1)


})
