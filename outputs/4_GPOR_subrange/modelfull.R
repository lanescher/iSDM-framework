code <- nimble::nimbleCode({


# Process Model

# 1177 cells
for (i in 1:nCell) {

    # log intensity
    XB0[i] <- inprod(B[1:nCovZ], Xz[i,1:nCovZ])
    log(lambda0[i]) <- XB0[i] + spat[i] # residual spatial effect

}


# ---------------------------------------

# Observation Model 1: PO, iNaturalist
# 1177 cells
for (j in 1:nW1) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD1[j] <- lambda0[Wcells1[j]] * alpha[1]

  # Observation model
  W1[j] ~ dpois(lambdaD1[j] * E1[j]) # Poisson

  # Effort 
  logit(eff1[j]) <- inprod(A1[1:nCovW1], Xw1[j,1:nCovW1])
  E1[j] <- eff1[j] * S1[j]

  # Prior for X imputation
  for (c in 1:nCovW1) {
    Xw1[j,c] ~ dnorm(0, 1)
  }
      
}

  

# Observation Model 2: PO, MD PO
# 1078 cells
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

  

# Observation Model 3: count, MDMBSS
# 198 observations, 1 median visits per site
for (j in 1:nY3) {

  # Make dataset-specific lambda (and N and Z if needed)
  lambdaD3[j] <- lambda0[Ycells3[j]] * alpha[3]

  # Observation model
  Y3[j] ~ dpois(lambdaD3[j] * p3[j]) # Poisson

  # Detection 
  log(p3[j]) <- inprod(C3[1:nCovY3], Xy3[j,1:nCovY3])

  # Prior for X imputation
  for (c in 1:nCovY3) {
    Xy3[j,c] ~ dnorm(0, 1)
  }
}

  

# Observation Model 4: count, NE Stream VES
# 23 observations, 1 median visits per site
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

# Observation priors, PO 2: MD PO
for (b in 1:nCovW2) {
  A2[b] ~ dnorm(0,1)
}

# Observation priors, count 3: MDMBSS
for (b in 1:nCovY3) {
  C3[b] ~ dnorm(0,1)
}

# Observation priors, count 4: NE Stream VES
for (b in 1:nCovY4) {
  C4[b] ~ dnorm(0,1)
}

tau <- 1
spat[1:nCell] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nCell], tau, zero_mean = 1)


})
