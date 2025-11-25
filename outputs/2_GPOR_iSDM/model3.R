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
# 49791 cells
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

  

# Observation Model 2: PO, MA PO, MD PO, ME PO, NY Atlas, VT PO, WV PO
# 49791 cells
for (j in 1:nW2) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD2[j] <- lambda0[Wcells2[j]] * alpha[2]

  # Observation model
  W2[j] ~ dpois(lambdaD2[j] * E2[j]) # Poisson

  # Effort 
  log(eff2[j]) <- inprod(A2[1:nCovW2], Xw2[j,1:nCovW2])
  E2[j] <- eff2[j] * S2[j]
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

  

# Observation Model 4: count, MA VES
# 179 observations, 1 median visits per site
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
# 1269 observations, 1 median visits per site
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
# 36 observations, 1 median visits per site
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
# 312 observations, 1 median visits per site
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

  

# Observation Model 8: count, Dodd Litterbag
# 133 observations, 9 median visits per site
for (j in 1:nY8) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD8[j] <- lambda0[Ycells8[j]] * alpha[8]

  # Observation model
  Y8[j] ~ dbinom(p8[j], ND8[j]) # N-mixture

  # Detection 
  logit(p8[j]) <- inprod(C8[1:nCovY8], Xy8[j,1:nCovY8])

  # Prior for X imputation
  for (c in 1:nCovY8) {
    Xy8[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 9: count, Dodd Terrestrial
# 132 observations, 1 median visits per site
for (j in 1:nY9) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD9[j] <- lambda0[Ycells9[j]] * alpha[9]

  # Observation model
  Y9[j] ~ dpois(lambdaD9[j] * p9[j]) # Poisson

  # Detection 
  log(p9[j]) <- C9[1] * Xy9[j,1]

  # Prior for X imputation
  Xy9[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 10: count, Hyde
# 929 observations, 9 median visits per site
for (j in 1:nY10) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD10[j] <- lambda0[Ycells10[j]] * alpha[10]

  # Observation model
  Y10[j] ~ dbinom(p10, ND10[j]) # N-mixture

  
}
# Detection 
logit(p10) <- C10[1] * 1
  
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

# Observation priors, PO 2: MA PO, MD PO, ME PO, NY Atlas, VT PO, WV PO
for (b in 1:nCovW2) {
  A2[b] ~ dnorm(0,1)
}

# Observation priors, DND 3: Dodd Aquatic
for (b in 1:nCovV3) {
  D3[b] ~ dnorm(0,1)
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

# Observation priors, count 8: Dodd Litterbag
for (b in 1:nCovY8) {
  C8[b] ~ dnorm(0,1)
}

# Observation priors, count 9: Dodd Terrestrial
for (b in 1:nCovY9) {
  C9[b] ~ dnorm(0,1)
}

# Observation priors, count 10: Hyde
for (b in 1) {
  C10[b] ~ dnorm(0,1)
}

tau <- 1
spat[1:nSpatCell] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nSpatCell], tau, zero_mean = 1)


})
