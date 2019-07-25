
# SigmaOut.rds contains the lower triangle of the variance covariance matrix, i.e., 136 columns with 16 species
Sigma <- readRDS(file="SigmaOut.rds")
str(Sigma)
dim(Sigma)

postSigMeans <- apply(Sigma,2,mean)
length(postSigMeans)

S <- diag(16)
S[lower.tri(S, diag=TRUE)] <- postSigMeans

# Convert S (half matrix) to full matrix
X <- diag(16)
X[lower.tri(X, diag=TRUE)] <- postSigMeans
X <- X + t(X) - diag(diag(X))


# Convert Sigma to Rho
Rho <- array(NA, dim=c(dim(Sigma)[1], dim(X)[1], dim(X)[1]))
dim(Rho)
for(samp in 1:dim(Sigma)[1]){
  # Grab mcmc sample from Sigma
  sigsamp <- Sigma[samp,]
  # Convert S (half matrix) to full matrix
  Stemp <- diag(16) 
  Stemp[lower.tri(Stemp, diag=TRUE)] <- sigsamp 
  Stemp <- Stemp + t(Stemp) - diag(diag(Stemp)) 
  # Calculate correlation matrix from Sigma
  Rho[samp, , ] <- cov2cor(Stemp)
}


dim(Rho)
RhoMean <- apply(Rho,c(3,2),mean)
# Make upper triangular matrix set to zero 
RhoMean[upper.tri(RhoMean)] <- NA
# remove diagonal
RhoMean[upper.tri(RhoMean,diag=TRUE)] <- NA
RMean <- as.numeric(t(RhoMean))
RMean <- RMean[!is.na(RMean)]

# Lower CI
RhoCIL <- apply(Rho,c(3,2),quantile, 0.025)
# Make upper triangular matrix set to zero 
RhoCIL[upper.tri(RhoCIL)] <- NA
# remove diagonal
RhoCIL[upper.tri(RhoCIL,diag=TRUE)] <- NA
RCIL <- as.numeric(t(RhoCIL))
RCIL <- RCIL[!is.na(RCIL)]

# Upper CI
RhoCIU <- apply(Rho,c(3,2),quantile, 0.975)
# Make upper triangular matrix set to zero 
RhoCIU[upper.tri(RhoCIU)] <- NA
# remove diagonal
RhoCIU[upper.tri(RhoCIU,diag=TRUE)] <- NA
RCIU <- as.numeric(t(RhoCIU))
RCIU <- RCIU[!is.na(RCIU)]

saveRDS(Rho, file="ResRho.rds")





