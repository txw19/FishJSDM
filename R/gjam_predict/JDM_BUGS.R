# rm(list=ls())
library(data.table)
library(lubridate)
library(car)
# library(jagsUI)
library(R2WinBUGS)
library(MCMCpack)
library(arm)
library(splitstackshape)
library(dplyr)
library(PerformanceAnalytics) 

dat <- fread('propagate_uncert.csv')
dat[, .N]
head(dat)
summary(dat)

# Remove observations with NA covariates
dat <- dat[!is.na(iws_urban)]
# dat <- dat[!is.na(hu12_dep_totaln_1990_mean)]
dat <- dat[!is.na(maxdepth)]

dim(dat)

# How many observations are very low for select response variables
sum(dat$median_secchi==0, na.rm=T)
sum(dat$median_tp <= 0.015, na.rm=T)
sum(dat$median_no2no3 <= 0, na.rm=T)

# Remove very low values for select response variables
dat <- dat[median_tp > 0.015 | is.na(median_tp)]
dat <- dat[median_secchi > 0.0 | is.na(median_secchi)]
dat <- dat[median_no2no3 > 0.0 | is.na(median_no2no3)]
# hist(log(dat$median_tp))

# log_transform response variables
dat[, log_tp:= log(median_tp)]
dat[,log_tn:=log(median_tn)]
dat[,log_chla:=log(median_chla)]
dat[, log_secchi:=log(median_secchi)]
dat[, log_no3:=log(median_no2no3)]

dim(dat)
head(dat)

# How many observations are NA
sum(is.na(dat$median_tp))
sum(is.na(dat$median_tn))
sum(is.na(dat$median_chla))
sum(is.na(dat$median_secchi))
sum(is.na(dat$median_no2no3))

# How many observations are observed
sum(!is.na(dat$median_tp))
sum(!is.na(dat$median_tn))
sum(!is.na(dat$median_chla))
sum(!is.na(dat$median_secchi))
sum(!is.na(dat$median_no2no3))


# Transform predictor variables
dat[,log_depth := scale(log(maxdepth))]
dat[,log_area := scale(log(lake_area_ha))]
dat[,log_wala := scale(log(wa_la))]
dat[, urban := scale(car::logit(iws_urban))]
dat[, ag := scale(car::logit(iws_ag))]
## dat[, pasture_ag := scale(car::logit(iws_pasture))]
## dat[, rowcrop_ag := scale(car::logit(iws_rowcrop))]
## dat[, forest := scale(car::logit(iws_total_forest))]
## dat[, decid_forest := scale(car::logit(iws_deciduous_forest))]
## dat[, evergreen_forest := scale(car::logit(iws_evergreen_forest))]
## dat[, mixed_forest := scale(car::logit(iws_mixed_forest))]
## dat[, woody_wetland := scale(car::logit(iws_woody_wetland))]
## dat[, emergent_wetland := scale(car::logit(iws_emergent_wetland))]
dat[, wetland := scale(car::logit(iws_total_wetland))]
dat[, road_density := scale(log(iws_roaddensity_density_mperha+0.01))]
dat[, stream_density := scale(log(iws_streamdensity_streams_density_mperha+0.01))]

summary(dat)
# cor(cbind(dat$log_depth, dat$log_area, dat$ag))
# cor(cbind(dat$log_depth, dat$log_area, dat$log_wala))

# Look at correlations among potential covariates
# chart.Correlation(dat[,71:86], histogram=TRUE, pch="+")
# chart.Correlation(dat[,82:86], histogram=TRUE, pch="+")
# Use total ag, total forest, total wetland based on correlations - not pasture, rowcrop, deciduous, evergreen, mixed
#                         emergen, woody
# Forest and Ag correlated (0.70) just use Ag

###### Create the out-of-sample test data set #######
# Index for splitting dat into 10 oos datasets (~10% of the total), will allow each lake will get predicted oos
ind <- rep(1:10,dim(dat)[1]/10)
ind <- c(ind, c(1,2,3))
length(ind)
# Randomize index
ind <- sample(ind)

dat$index <- ind
table(dat$index)
summary(dat)
# max(table(dat$index,dat$lagoslakeid))
# write out final dataset 
write.csv(dat, "dat.csv", row.names = F)
dim(dat)

########### Begin loop through datasets for 10-fold cross-validation
for(i in 1:10){

# Select out dataset for model fitting (datBugs)
datBugs <- dat[index != i,]
dim(datBugs)
# summary(datBugs)
# table(datBugs$index)  
  


#################################################################
########## BUGS CODE ############################################
#################################################################

# Define the model in the BUGS language and write a text file
sink("jdm.txt")
cat("
    model {
    
for(i in 1:n){
  y[i,1:K] ~ dmnorm(mu[i,], Tau.B[,])
    for(j in 1:5){
    mu[i,j] <- alpha[j] + beta[j,1]*x[i,1] + beta[j,2]*x[i,2] + beta[j,3]*x[i,3] + beta[j,4]*x[i,4] + beta[j,5]*x[i,5]
              + beta[j,6]*x[i,6] + beta[j,7]*x[i,7] + beta[j,8]*x[i,8]
    }
  }


# Model variance-covariance
Tau.B[1:K,1:K] ~ dwish(W[,], df)
df <- K+1
Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
for (k in 1:K){
  for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
  }
  sigma.B[k] <- sqrt(Sigma.B[k,k])
}

## Priors for regression coefficients
for(j in 1:K){
    alpha[j] ~ dnorm(0, 0.0001)
   for(k in 1:nbeta) {
     beta[j,k] ~ dnorm(0, 0.0001)
   }
}
    

} # end model
    ",fill = TRUE)
sink()


head(datBugs)

# Number of response variables
K <- 5



# Create identity matrix for Wishart dist'n
W <- diag(K)

head(datBugs)

# Create response matrix
y_mat <- cbind(datBugs$log_tp, datBugs$log_tn, datBugs$log_chla, datBugs$log_secchi, datBugs$log_no3)
dim(y_mat)


# Predictor matrix
# x <- cbind(datBugs$log_depth, datBugs$log_area)
x <- as.matrix(datBugs[,71:78])
# Number of predictors
J <- dim(x)[2]

# chart.Correlation(dat[,71:79], histogram=TRUE, pch="+")

# Load data
data <- list(y = y_mat, n = dim(datBugs)[1], x=x,K=K,W=W,nbeta=dim(x)[2] )


# Initial values
# inits <- function (){
#   list (alpha = rnorm(K),Tau.B=rwish(K+1,diag(K)) )
# }
inits <- function (){
  list (alpha = rnorm(K),
        beta = matrix(rnorm(J*K),nrow=K,ncol=J),
        Tau.B=rwish(K+1,diag(K)) )
}

# Parameters monitored
parameters <- c("alpha","beta","Sigma.B","rho.B")

# ni = 15000
# nb = 5000
# MCMC settings
ni <- 15000
nt <- 1
nb <- 5000
nc <- 3

bugs.dir <- "C:/Program Files/WinBUGS14/"
 
start.time = Sys.time()         # Set timer 
out <- bugs(data = data, inits = inits, parameters.to.save = parameters, 
            model.file = "jdm.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb,debug = F, bugs.directory=bugs.dir, DIC=FALSE)
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time

# Two predictors = 1.7 run time

# Summarize posteriors
# print(out, dig = 3)

# arm::traceplot(out)

# Export model fit summary for each model fit
BugsOut <- out$summary
write.csv(BugsOut, paste("BUGSModelSummary",i,".csv", sep=''), row.names = T)
# Export MCMC samples for each model fit
mcmcOut <- out$sims.list
saveRDS(mcmcOut, file=paste("JDMmcmc_out",i,".rds",sep=''))

} # end for loop


