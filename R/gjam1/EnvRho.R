# rm(list=ls())
require(data.table)
require(mvtnorm)

### Caclulate, from a JDM, the between species correlation due to their shared
### environmental responses. Following the approach of Pollock et al. (2014).

# Predict Y1 - Y5 (log_tp, log_tn, log_chla, log_secchi, log_no3)

# Read in final dat
datFinal <- fread("dat.csv")

# Read in RDS file model fit
out1 <- readRDS(file="JDMmcmc_out1.rds")

# Number of response variables
numY <- dim(out1$alpha)[2]

# Grab X's from datFinal (needs to be !=1 form .rds out1)
datOOS <- datFinal[index != 1,] 
####### !!!!!!!! Need to specify predictors
# newX <- cbind(datOOS$log_depth, datOOS$log_area)
X <- as.matrix(datOOS[,71:78])
# Add 1's to design matrix for matrix multiplication (not needed for this calc, but added for adapting Pollock code)
X <- as.matrix(cbind(c(rep(1,dim(X)[1])),X))

# Number of covariates
K <- dim(X[,-1])[2]

# number of sites 
n.sites <- nrow(X)

# number of species in Occurance matrix
n.species <- numY

# Beta coefficients
Beta <- out1$beta


# Number of mcmc samples 
n.sims <- dim(Beta)[1]

# Variance covariance matrix (named Sigma2 for adpating Pollock code)
Sigma2 <- out1$Sigma.B

############# Calculate the correlation due to the environment
EnvRho <- apply(Beta, 1, 
                function(x) {
                  matrix(rowSums(apply(cbind(x[, -1]), 2, function(y) outer(y, y))), n.species)
                }
)
dim(EnvRho) <- rev(dim(Sigma2))
EnvRho <- aperm(EnvRho, c(3, 2, 1))

COVX <- cov(X)

for(sims in seq_len(n.sims)) {
  for(species in seq_len(n.species)) {
    for(species.prime in seq_len(n.species)) {
      EnvRho[sims, species, species.prime] <- (EnvRho[sims, species, species.prime] +
                                                 sum(sapply(seq_len(K)[-1], function(k) {
                                                   sum(Beta[sims, species, k] * 
                                                         Beta[sims, species.prime, seq_len(K)[-c(1, k)]] * 
                                                         COVX[k, seq_len(K)[-c(1, k)]]
                                                   )}
                                                 )))
    }
  }
}

EnvRho <- apply(EnvRho, 1, cov2cor)
dim(EnvRho) <- rev(dim(Sigma2))
EnvRho <- aperm(EnvRho, c(3, 2, 1))

dim(EnvRho)

# Posterior mean EnvRho and 95% credile intervals
apply(EnvRho,c(3,2),mean)
apply(EnvRho,c(3,2),quantile, 0.025)
apply(EnvRho,c(3,2),quantile, 0.975)

# Save posterior samples for EnvRho
saveRDS(EnvRho, file="EnvRho.rds")



