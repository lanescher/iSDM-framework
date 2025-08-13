code <- nimbleCode({


# Process Model

# 25616 cells
for (i in 1:nCell) {

    # log intensity
    XB0[i] <- inprod(B[1:nCovZ], Xz[i,1:nCovZ])
    log(lambda0[i]) <- XB0[i] + spat[i] # residual spatial effect

}


# ---------------------------------------

# Observation Model 1: PO, iNaturalist
# 16940 cells
for (j in 1:nW1) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD1[j] <- lambda0[Wcells1[j]] * alpha[1]

  # Observation model
  W1[j] ~ dpois(lambdaD1[j] * E1[j]) # Poisson

  # Effort 
  logit(eff1[j]) <- A1[1] * Xw1[j,1]
  E1[j] <- eff1[j] * S1[j]

  # Prior for X imputation
  Xw1[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 2: PO, Museum
# 16940 cells
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

  

# Observation Model 3: count, USGS VES - MS
# 91 observations, 1 median visits per site
for (j in 1:nY3) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD3[j] <- lambda0[Ycells3[j]] * alpha[3]

  # Observation model
  Y3[j] ~ dpois(lambdaD3[j] * p3) # Poisson

  
}
# Detection 
log(p3) <- 0
  

# Observation Model 4: count, Odonnell
# 1078 observations, 27 median visits per site
for (j in 1:nY4) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD4[j] <- lambda0[Ycells4[j]] * alpha[4]

  # Observation model
  Y4[j] ~ dbinom(p4, ND4[j]) # N-mixture

  
}
# Detection 
logit(p4) <- C4[1] * 1
  
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

# Observation priors, PO 1: iNaturalist
for (b in 1:nCovW1) {
  A1[b] ~ dnorm(0,1)
}

# Observation priors, PO 2: Museum
for (b in 1:nCovW2) {
  A2[b] ~ dnorm(0,1)
}


# Observation priors, count 4: Odonnell
for (b in 1) {
  C4[b] ~ dnorm(0,1)
}

tau <- 1
spat[1:nCell] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nCell], tau, zero_mean = 1)


})
