code <- nimbleCode({


# Process Model

# 76016 cells
for (i in 1:nCell) {

    # log intensity
    XB0[i] <- inprod(B[1:nCovZ], Xz[i,1:nCovZ])
    log(lambda0[i]) <- XB0[i] + spat[spatCells[i]] # residual spatial effect

}


# ---------------------------------------

# Observation Model 1: PO, iNaturalist
# 57842 cells
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
# 57842 cells
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

  

# Observation Model 3: PO, NY Atlas, WV PO, VT PO, MA PO
# 57842 cells
for (j in 1:nW3) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD3[j] <- lambda0[Wcells3[j]] * alpha[3]

  # Observation model
  W3[j] ~ dpois(lambdaD3[j] * E3[j]) # Poisson

  # Effort 
  log(eff3[j]) <- inprod(A3[1:nCovW3], Xw3[j,1:nCovW3])
  E3[j] <- eff3[j] * S3[j]
}

  

# Observation Model 4: DND, Dodd Aquatic
# 158 observations, 1 median visits per site
for (j in 1:nV4) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD4[j] <- lambda0[Vcells4[j]] * alpha[4]

  # Observation model
  V4[j] ~ dbern(1-exp(-lambdaD4[j] * p4[j])) # Bernoulli

  # Detection 
  log(p4[j]) <- D4[1] * Xv4[j,1]

  # Prior for X imputation
  Xv4[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 5: count, Hyde
# 938 observations, 9 median visits per site
for (j in 1:nY5) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD5[j] <- lambda0[Ycells5[j]] * alpha[5]

  # Observation model
  Y5[j] ~ dbinom(p5, ND5[j]) # N-mixture

  
}
# Detection 
logit(p5) <- C5[1] * 1
  

# Observation Model 6: count, Dodd Terrestrial
# 132 observations, 1 median visits per site
for (j in 1:nY6) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD6[j] <- lambda0[Ycells6[j]] * alpha[6]

  # Observation model
  Y6[j] ~ dpois(lambdaD6[j] * p6[j]) # Poisson

  # Detection 
  log(p6[j]) <- C6[1] * Xy6[j,1]

  # Prior for X imputation
  Xy6[j, 1] ~ dnorm(0, 1)
}

  

# Observation Model 7: count, Dodd Litterbag
# 133 observations, 9 median visits per site
for (j in 1:nY7) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD7[j] <- lambda0[Ycells7[j]] * alpha[7]

  # Observation model
  Y7[j] ~ dbinom(p7[j], ND7[j]) # N-mixture

  # Detection 
  logit(p7[j]) <- inprod(C7[1:nCovY7], Xy7[j,1:nCovY7])

  # Prior for X imputation
  for (c in 1:nCovY7) {
    Xy7[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 8: count, NE Stream VES
# 256 observations, 1 median visits per site
for (j in 1:nY8) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD8[j] <- lambda0[Ycells8[j]] * alpha[8]

  # Observation model
  Y8[j] ~ dpois(lambdaD8[j] * p8[j]) # Poisson

  # Detection 
  log(p8[j]) <- inprod(C8[1:nCovY8], Xy8[j,1:nCovY8])

  # Prior for X imputation
  for (c in 1:nCovY8) {
    Xy8[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 9: count, MA VES
# 179 observations, 1 median visits per site
for (j in 1:nY9) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD9[j] <- lambda0[Ycells9[j]] * alpha[9]

  # Observation model
  Y9[j] ~ dpois(lambdaD9[j] * p9[j]) # Poisson

  # Detection 
  log(p9[j]) <- inprod(C9[1:nCovY9], Xy9[j,1:nCovY9])

  # Prior for X imputation
  for (c in 1:nCovY9) {
    Xy9[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 10: count, MDMBSS
# 1378 observations, 1 median visits per site
for (j in 1:nY10) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD10[j] <- lambda0[Ycells10[j]] * alpha[10]

  # Observation model
  Y10[j] ~ dpois(lambdaD10[j] * p10[j]) # Poisson

  # Detection 
  log(p10[j]) <- inprod(C10[1:nCovY10], Xy10[j,1:nCovY10])

  # Prior for X imputation
  for (c in 1:nCovY10) {
    Xy10[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 11: count, SE Stream VES
# 16 observations, 1 median visits per site
for (j in 1:nY11) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD11[j] <- lambda0[Ycells11[j]] * alpha[11]

  # Observation model
  Y11[j] ~ dpois(lambdaD11[j] * p11[j]) # Poisson

  # Detection 
  log(p11[j]) <- inprod(C11[1:nCovY11], Xy11[j,1:nCovY11])

  # Prior for X imputation
  for (c in 1:nCovY11) {
    Xy11[j,c] ~ dnorm(0, 1)
  }
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

# Observation priors, PO 3: NY Atlas, WV PO, VT PO, MA PO
for (b in 1:nCovW3) {
  A3[b] ~ dnorm(0,1)
}

# Observation priors, DND 4: Dodd Aquatic
for (b in 1:nCovV4) {
  D4[b] ~ dnorm(0,1)
}

# Observation priors, count 5: Hyde
for (b in 1) {
  C5[b] ~ dnorm(0,1)
}

# Observation priors, count 6: Dodd Terrestrial
for (b in 1:nCovY6) {
  C6[b] ~ dnorm(0,1)
}

# Observation priors, count 7: Dodd Litterbag
for (b in 1:nCovY7) {
  C7[b] ~ dnorm(0,1)
}

# Observation priors, count 8: NE Stream VES
for (b in 1:nCovY8) {
  C8[b] ~ dnorm(0,1)
}

# Observation priors, count 9: MA VES
for (b in 1:nCovY9) {
  C9[b] ~ dnorm(0,1)
}

# Observation priors, count 10: MDMBSS
for (b in 1:nCovY10) {
  C10[b] ~ dnorm(0,1)
}

# Observation priors, count 11: SE Stream VES
for (b in 1:nCovY11) {
  C11[b] ~ dnorm(0,1)
}

tau ~ dgamma(0.01,0.01)
spat[1:nSpatCell] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nSpatCell], tau, zero_mean = 1)


})
