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
# 13804 cells
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

  

# Observation Model 2: PO, FS NRIS
# 13804 cells
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

  

# Observation Model 3: PO, ORBIC
# 13804 cells
for (j in 1:nW3) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD3[j] <- lambda0[Wcells3[j]] * alpha[3]

  # Observation model
  W3[j] ~ dpois(lambdaD3[j] * E3[j]) # Poisson

  # Effort 
  log(E3[j]) <- A3[1] * Xw3[j,1]

  # Prior for X imputation
  Xw3[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 4: PO, WDWF
# 13804 cells
for (j in 1:nW4) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD4[j] <- lambda0[Wcells4[j]] * alpha[4]

  # Observation model
  W4[j] ~ dpois(lambdaD4[j] * E4[j]) # Poisson

  # Effort 
  log(E4[j]) <- A4[1] * Xw4[j,1]

  # Prior for X imputation
  Xw4[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 5: PO, NCCN
# 13804 cells
for (j in 1:nW5) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD5[j] <- lambda0[Wcells5[j]] * alpha[5]

  # Observation model
  W5[j] ~ dpois(lambdaD5[j] * E5[j]) # Poisson

  # Effort 
  log(E5[j]) <- A5[1] * Xw5[j,1]

  # Prior for X imputation
  Xw5[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 6: PO, Museum
# 13804 cells
for (j in 1:nW6) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD6[j] <- lambda0[Wcells6[j]] * alpha[6]

  # Observation model
  W6[j] ~ dpois(lambdaD6[j] * E6[j]) # Poisson

  # Effort 
  log(E6[j]) <- A6[1] * Xw6[j,1]

  # Prior for X imputation
  Xw6[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 7: DND, USGS RACA assessment - historic
# 48 observations, 1 median visits per site
for (j in 1:nV7) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD7[j] <- lambda0[Vcells7[j]] * alpha[7]

  # Observation model
  V7[j] ~ dbern(1-exp(-lambdaD7[j] * p7[j])) # Bernoulli

  # Detection 
  log(p7[j]) <- inprod(D7[1:nCovV7], Xv7[j,1:nCovV7])

  # Prior for X imputation
  for (c in 1:nCovV7) {
    Xv7[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 8: count, USGS Count
# 26 observations, 2 median visits per site
for (j in 1:nY8) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD8[j] <- lambda0[Ycells8[j]] * alpha[8]

  # Observation model
  Y8[j] ~ dpois(lambdaD8[j] * p8[j]) # Poisson

  # Detection 
  log(p8[j]) <- C8[1] * Xy8[j,1]

  # Prior for X imputation
  Xy8[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 9: count, ARMI VES
# 334 observations, 11 median visits per site
for (j in 1:nY9) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD9[j] <- lambda0[Ycells9[j]] * alpha[9]

  # Observation model
  Y9[j] ~ dbinom(p9[j], ND9[j]) # N-mixture

  # Detection 
  logit(p9[j]) <- inprod(C9[1:nCovY9], Xy9[j,1:nCovY9])

  # Prior for X imputation
  for (c in 1:nCovY9) {
    Xy9[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 10: count, USGS RACA assessment - current
# 42 observations, 1 median visits per site
for (j in 1:nY10) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD10[j] <- lambda0[Ycells10[j]] * alpha[10]

  # Observation model
  Y10[j] ~ dpois(lambdaD10[j] * p10[j]) # Poisson

  # Detection 
  log(p10[j]) <- C10[1] * Xy10[j,1]

  # Prior for X imputation
  Xy10[j, 1] ~ dnorm(0, 1)
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

# Observation priors, PO 1: BLM Surveys
for (b in 1:nCovW1) {
  A1[b] ~ dnorm(0,1)
}

# Observation priors, PO 2: FS NRIS
for (b in 1:nCovW2) {
  A2[b] ~ dnorm(0,1)
}

# Observation priors, PO 3: ORBIC
for (b in 1:nCovW3) {
  A3[b] ~ dnorm(0,1)
}

# Observation priors, PO 4: WDWF
for (b in 1:nCovW4) {
  A4[b] ~ dnorm(0,1)
}

# Observation priors, PO 5: NCCN
for (b in 1:nCovW5) {
  A5[b] ~ dnorm(0,1)
}

# Observation priors, PO 6: Museum
for (b in 1:nCovW6) {
  A6[b] ~ dnorm(0,1)
}

# Observation priors, DND 7: USGS RACA assessment - historic
for (b in 1:nCovV7) {
  D7[b] ~ dnorm(0,1)
}

# Observation priors, count 8: USGS Count
for (b in 1:nCovY8) {
  C8[b] ~ dnorm(0,1)
}

# Observation priors, count 9: ARMI VES
for (b in 1:nCovY9) {
  C9[b] ~ dnorm(0,1)
}

# Observation priors, count 10: USGS RACA assessment - current
for (b in 1:nCovY10) {
  C10[b] ~ dnorm(0,1)
}

tau <- 1
spat[1:nCell] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nCell], tau, zero_mean = 1)


})
