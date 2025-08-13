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
# 25616 cells
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
# 25616 cells
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

  

# Observation Model 3: DND, Dodd Aquatic
# 158 observations, 1 median visits per site
for (j in 1:nV3) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD3[j] <- lambda0[Vcells3[j]] * alpha[3]

  # Observation model
  V3[j] ~ dbern(1-exp(-lambdaD3[j] * p3[j])) # Bernoulli

  # Detection 
  log(p3[j]) <- D3[1] * Xv3[j,1]

  # Prior for X imputation
  Xv3[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 4: count, USGS VES - MS
# 91 observations, 1 median visits per site
for (j in 1:nY4) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD4[j] <- lambda0[Ycells4[j]] * alpha[4]

  # Observation model
  Y4[j] ~ dpois(lambdaD4[j] * p4) # Poisson

  
}
# Detection 
log(p4) <- 0
  

# Observation Model 5: count, Odonnell
# 1078 observations, 27 median visits per site
for (j in 1:nY5) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD5[j] <- lambda0[Ycells5[j]] * alpha[5]

  # Observation model
  Y5[j] ~ dbinom(p5, ND5[j]) # N-mixture

  
}
# Detection 
logit(p5) <- C5[1] * 1
  

# Observation Model 6: count, Hyde
# 929 observations, 9 median visits per site
for (j in 1:nY6) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD6[j] <- lambda0[Ycells6[j]] * alpha[6]

  # Observation model
  Y6[j] ~ dbinom(p6, ND6[j]) # N-mixture

  
}
# Detection 
logit(p6) <- C6[1] * 1
  

# Observation Model 7: count, Bailey
# 384 observations, 10 median visits per site
for (j in 1:nY7) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD7[j] <- lambda0[Ycells7[j]] * alpha[7]

  # Observation model
  Y7[j] ~ dbinom(p7, ND7[j]) # N-mixture

  
}
# Detection 
logit(p7) <- C7[1] * 1
  

# Observation Model 8: count, Dodd Terrestrial
# 132 observations, 1 median visits per site
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

# Observation priors, DND 3: Dodd Aquatic
for (b in 1:nCovV3) {
  D3[b] ~ dnorm(0,1)
}


# Observation priors, count 5: Odonnell
for (b in 1) {
  C5[b] ~ dnorm(0,1)
}

# Observation priors, count 6: Hyde
for (b in 1) {
  C6[b] ~ dnorm(0,1)
}

# Observation priors, count 7: Bailey
for (b in 1) {
  C7[b] ~ dnorm(0,1)
}

# Observation priors, count 8: Dodd Terrestrial
for (b in 1:nCovY8) {
  C8[b] ~ dnorm(0,1)
}

tau <- 1
spat[1:nCell] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nCell], tau, zero_mean = 1)


})
