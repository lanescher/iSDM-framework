code <- nimbleCode({


# Process Model

# 66011 cells
for (i in 1:nCell) {

    # log intensity
    XB0[i] <- inprod(B[1:nCovZ], Xz[i,1:nCovZ])
    log(lambda0[i]) <- XB0[i] + spat[spatCells[i]] # residual spatial effect

}


# ---------------------------------------

# Observation Model 1: PO, iNaturalist
# 29337 cells
for (j in 1:nW1) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD1[j] <- lambda0[Wcells1[j]] * alpha[1]

  # Observation model
  W1[j] ~ dpois(lambdaD1[j] * E1[j]) # Poisson

  # Effort 
  logit(eff1[j]) <- inprod(A1[1:nCovW1], Xw1[j,1:nCovW1])
  E1[j] <- eff1[j] * S1[j]

  # Prior for X imputation
  Xw1[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 2: PO, Museum
# 29337 cells
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

  

# Observation Model 3: PO, ME PO, NY Atlas, VT PO, WV PO
# 29337 cells
for (j in 1:nW3) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD3[j] <- lambda0[Wcells3[j]] * alpha[3]

  # Observation model
  W3[j] ~ dpois(lambdaD3[j] * E3[j]) # Poisson

  # Effort 
  log(eff3[j]) <- inprod(A3[1:nCovW3], Xw3[j,1:nCovW3])
  E3[j] <- eff3[j] * S3[j]
}

  

# Observation Model 4: count, MA VES
# 8 observations, 1 median visits per site
for (j in 1:nY4) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD4[j] <- lambda0[Ycells4[j]] * alpha[4]

  # Observation model
  Y4[j] ~ dpois(lambdaD4[j] * p4[j]) # Poisson

  # Detection 
  log(p4[j]) <- inprod(C4[1:nCovY4], Xy4[j,1:nCovY4])

  # Prior for X imputation
  for (c in 1:nCovY4) {
    Xy4[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 5: count, MDMBSS
# 63 observations, 1 median visits per site
for (j in 1:nY5) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD5[j] <- lambda0[Ycells5[j]] * alpha[5]

  # Observation model
  Y5[j] ~ dpois(lambdaD5[j] * p5[j]) # Poisson

  # Detection 
  log(p5[j]) <- inprod(C5[1:nCovY5], Xy5[j,1:nCovY5])

  # Prior for X imputation
  for (c in 1:nCovY5) {
    Xy5[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 6: count, SE Stream VES
# 22 observations, 1 median visits per site
for (j in 1:nY6) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD6[j] <- lambda0[Ycells6[j]] * alpha[6]

  # Observation model
  Y6[j] ~ dpois(lambdaD6[j] * p6[j]) # Poisson

  # Detection 
  log(p6[j]) <- inprod(C6[1:nCovY6], Xy6[j,1:nCovY6])

  # Prior for X imputation
  for (c in 1:nCovY6) {
    Xy6[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 7: count, NE Stream VES
# 140 observations, 1 median visits per site
for (j in 1:nY7) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD7[j] <- lambda0[Ycells7[j]] * alpha[7]

  # Observation model
  Y7[j] ~ dpois(lambdaD7[j] * p7[j]) # Poisson

  # Detection 
  log(p7[j]) <- inprod(C7[1:nCovY7], Xy7[j,1:nCovY7])

  # Prior for X imputation
  for (c in 1:nCovY7) {
    Xy7[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 8: count, Hyde
# 25 observations, 8 median visits per site
for (j in 1:nY8) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD8[j] <- lambda0[Ycells8[j]] * alpha[8]

  # Observation model
  Y8[j] ~ dbinom(p8, ND8[j]) # N-mixture

  
}
# Detection 
logit(p8) <- C8[1] * 1
  
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

# Observation priors, PO 3: ME PO, NY Atlas, VT PO, WV PO
for (b in 1:nCovW3) {
  A3[b] ~ dnorm(0,1)
}

# Observation priors, count 4: MA VES
for (b in 1:nCovY4) {
  C4[b] ~ dnorm(0,1)
}

# Observation priors, count 5: MDMBSS
for (b in 1:nCovY5) {
  C5[b] ~ dnorm(0,1)
}

# Observation priors, count 6: SE Stream VES
for (b in 1:nCovY6) {
  C6[b] ~ dnorm(0,1)
}

# Observation priors, count 7: NE Stream VES
for (b in 1:nCovY7) {
  C7[b] ~ dnorm(0,1)
}

# Observation priors, count 8: Hyde
for (b in 1) {
  C8[b] ~ dnorm(0,1)
}

tau <- 1
spat[1:nSpatCell] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nSpatCell], tau, zero_mean = 1)


})
