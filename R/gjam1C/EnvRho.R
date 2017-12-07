# rm(list=ls())
require(data.table)
require(mvtnorm)
require(abind)

### Caclulate, from a JDM, the between species correlation due to their shared
### environmental responses. Following the approach of Pollock et al. (2014).


# Read in final dat
datFinal <- fread("Finaldat.csv")

# Read in RDS file model fit
out1 <- readRDS(file="gjamOUT1.rds")

# Number of response variables
numY <- 16

####### !!!!!!!! Need to specify predictors

X <- as.matrix(scale(datFinal[,20:24]))
# Add 1's to design matrix for matrix multiplication (not needed for this calc, but added for adapting Pollock code)
X <- as.matrix(cbind(c(rep(1,dim(X)[1])),X))
dim(X)
head(X)


# Number of covariates
K <- dim(X[,-1])[2]

# number of sites 
n.sites <- nrow(X)

# number of species in Occurance matrix
n.species <- numY

# Beta coefficients dim = [nsim, 96], need array of [nsim, nspecies, ncovariates]
Beta <- out1$chains$bgibbs
dim(Beta)

# Remove intercepts from Beta

int.names <- c("black.bullhead_intercept","black.crappie_intercept","bowfin.(dogfish)_intercept","brown.bullhead_intercept",
  "common.carp_intercept","golden.shiner_intercept","largemouth.bass_intercept","northern.pike_intercept",
  "rock.bass_intercept","smallmouth.bass_intercept","tullibee.(cisco)_intercept","walleye_intercept",
  "white.sucker_intercept","yellow.bullhead_intercept","yellow.perch_intercept","sunfish_intercept")

Beta2 <- Beta[, !colnames(Beta) %in% int.names]

dim(Beta2)
head(Beta2)

# apply(Beta2,2,mean)

# NUmber of MCMC samples to use for calculating EnvRho
nsim <- 50000
# Chain length from analysis
chainLength <- dim(Beta2)[1]
# # Select thinned steps in chain for posterior predictions to ensure we take values from 
# # length of posterior
ID = seq( 1 , chainLength , floor(chainLength/nsim) )
length(ID)
# Grab nsim samples
Beta2sub <- Beta2[ID,]
dim(Beta2sub)
head(Beta2sub)

# We need matrices, one for each covariate, with rows=num MCMC samples, columns = number of species
# and a list for each covariates. 
const <- 5
lst1 <- list()
for(i in 1:5){
    lst1[[i]] <- Beta2sub[,c(i,5+i,10+i,15+i,20+i,25+i,30+i,35+i,40+i,45+i,50+i,55+i,60+i,65+i,70+i,75+i)]
}
str(lst1)
head(lst1[[1]])

BetaA <- array(as.numeric(unlist(lst1)), dim=c(length(ID),numY,K))

dim(BetaA)

# Rename for code below
Beta <- BetaA


# Number of mcmc samples 
n.sims <- dim(Beta)[1]

# Variance covariance matrix (named Sigma2 for adpating Pollock code), just need this for dimensions

Sigma2 <- array(NA, c(length(ID), numY, numY))
dim(Sigma2)

############# Calculate the correlation due to the environment
EnvRho <- apply(Beta, 1, 
                function(x) {
                  matrix(rowSums(apply(cbind(x[, -1]), 2, function(y) outer(y, y))), n.species)
                }
)
dim(EnvRho) <- rev(dim(Sigma2))
EnvRho <- aperm(EnvRho, c(3, 2, 1))
dim(EnvRho)
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



