# Make projections
data$Xz1 <- data$Xz
data$Xz1$sqrtarea_small <- data$Xz1$sqrtarea_small + rnorm(length(data$Xz1$sqrtarea_small), 1, 1)

data$Xz2 <- data$Xz
data$Xz2$sqrtarea_small <- data$Xz2$sqrtarea_small + rnorm(length(data$Xz2$sqrtarea_small), 2, 1)

proj.names <- c("proj1", "proj2")


# End script