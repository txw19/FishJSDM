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
numY <- 17

####### !!!!!!!! Need to specify predictors
# area.hectares + max_depth_m + mean.gdd + proportion_disturbed_ws + littoral.zone + secchi.m

X <- as.matrix(scale(datFinal[,2:7]))
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
###! Need to read in and create ydat from 1_gjam_all_spp.R to run this


fishP= read.csv("finalDat.csv")

ydat <- fishP[,8:24]
ydat <- as.data.frame(ydat)
spp.names <- colnames(ydat)
int.names <- vector()
for(i in 1:dim(ydat)[2]){
  int.names[i] <- paste(spp.names[i],"_intercept",sep='')
}
Beta2 <- Beta[, !colnames(Beta) %in% int.names]
dim(Beta2)


dim(Beta2)
head(Beta2)


# apply(Beta2,2,mean)

# NUmber of MCMC samples to use for calculating EnvRho
nsim <- 5000
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
const <- 6 # This is the number of covariates, we are stepping through and creating, for instance, a matrix
# of the mcmc samples (rows) for the effects of stream width for each species (column)
lst1 <- list()
# for(i in 1:6){
#     lst1[[i]] <- Beta2sub[,c(i,6+i,11+i,16+i,21+i,26+i,31+i,36+i,41+i,46+i,51+i,56+i,61+i,66+i,71+i,76+i,81+i,86+i,91+i,96+i)]
# }
# str(lst1)
# head(lst1[[1]])

spp.names <- colnames(ydat)
#right now output has coldwaterspp
spp.names[17]="coldwaterspp"

param.names <- vector()
for(i in 1:dim(ydat)[2]){
  param.names[i] <- paste(spp.names[i],"_area.hectares",sep='')
}
lst1[[1]] <- Beta2sub[, colnames(Beta2sub) %in% param.names]
head(lst1[[1]])

param.names <- vector()
for(i in 1:dim(ydat)[2]){
  param.names[i] <- paste(spp.names[i],"_maxdepthm",sep='')
}
lst1[[2]] <- Beta2sub[, colnames(Beta2sub) %in% param.names]
head(lst1[[2]])

param.names <- vector()
for(i in 1:dim(ydat)[2]){
  param.names[i] <- paste(spp.names[i],"_mean.gdd",sep='')
}
lst1[[3]] <- Beta2sub[, colnames(Beta2sub) %in% param.names]
head(lst1[[3]])

param.names <- vector()
for(i in 1:dim(ydat)[2]){
  param.names[i] <- paste(spp.names[i],"_proportiondisturbedws",sep='')
}
lst1[[4]] <- Beta2sub[, colnames(Beta2sub) %in% param.names]
head(lst1[[4]])

param.names <- vector()
for(i in 1:dim(ydat)[2]){
  param.names[i] <- paste(spp.names[i],"_littoral.zone",sep='')
}
lst1[[5]] <- Beta2sub[, colnames(Beta2sub) %in% param.names]
head(lst1[[5]])

param.names <- vector()
for(i in 1:dim(ydat)[2]){
  param.names[i] <- paste(spp.names[i],"_secchi.m",sep='')
}
lst1[[6]] <- Beta2sub[, colnames(Beta2sub) %in% param.names]
head(lst1[[6]])
dim(lst1[[6]])

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



