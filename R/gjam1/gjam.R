# rm(list=ls())
# library(gjam)
library(data.table)


library('Rcpp')
library('RcppArmadillo')
# placed cppFns.cpp in working directory and it worked
sourceCpp('cppFns.cpp')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjam.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamPlot.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamSimData.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamTrimY.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamCensorY.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamDeZero.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamReZero.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamPredict.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamSpec2Trait.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamIIE.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamIIEplot.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamOrdination.R')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamPoints2Grid.R')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamPriorTemplate.r')
source('C:/Users/txw19/Documents/Manuscripts/Joint_Fish_Modeling/gjam/R/gjamHfunctions.r')
library(MASS)


# Read in lake covariate and fish data
lakeP=fread("lake_predictors_for_joint_distribution.csv")
fishP=fread("cpue_JDM.csv")
dim(lakeP) #1343 x 8
dim(fishP) #1343 x 17

lakeP
fishP

# Look at correlations among covariates
# cor(cbind(lakeP$area.hectares, lakeP$max_depth_m, lakeP$Secchi.lake.mean,lakeP$mean.gdd, lakeP$alkalinity))

# transform covariates (covariates are standardized by mean and variance, then transformed back to original scales in output in gjam)
lakeP[,log_area := log(area.hectares)]
lakeP[,log_depth := log(max_depth_m)]
# lakeP[,z_secchi := scale(Secchi.lake.mean)]
# lakeP[,z_gdd := scale(mean.gdd)]
# lakeP[, z_alk := scale(alkalinity)]

# Merge cpe and predictors
setkey(lakeP, DOW)
setkey(fishP, DOW)
dat <- merge(fishP, lakeP)
dat[,.N]
dim(dat)

write.csv(dat, "Finaldat.csv", row.names = F)

dim(dat)

xdat <- dat[,20:24]
ydat <- sqrt(dat[,2:17])

xdat <- as.data.frame(xdat)
ydat <- as.data.frame(ydat)


# ml <-list(PREDICTX = F, ng=1000,burnin=500,typeNames=rep("CA",16)) 
# ml$FULL=T
# 
# jdmtest = gjam(~ log_area + log_depth + Secchi.lake.mean + mean.gdd + alkalinity, xdata=xdat, ydata=ydat,
#             modelList=ml)


# Some gjam options for modelList:
# holdoutN = 0, number of observations to hold out for out-of-sample prediction.
# holdoutIndex =  numeric(0), numeric vector of observations (row numbers) to holdout for out-of-sample prediction
# ,holdoutN=200
# -typeNames = "CA" means continuous abundance (meaning greater than 0).
# modelList$FULL = T
# start.time = Sys.time()  # Start timer 

ml <-list(PREDICTX = F, ng=50000,burnin=30000,typeNames=rep("CA",16)) 
ml$FULL=T

jdm1 = gjam(~ log_area + log_depth + Secchi.lake.mean + mean.gdd + alkalinity, xdata=xdat, ydata=ydat,
          modelList=ml)

# end.time = Sys.time()
# elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
# cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

# jdm1 <- readRDS(file="gjamOUT1.rds")
str(jdm1)
summary(jdm1)
names(jdm1)

# in-sample posterior predictions
fullpostY <- jdm1$chains$ygibbs
dim(fullpostY)
yMedian_gjam <- apply(fullpostY,2,median)
length(yMedian_gjam)

yMean_gjam <- apply(fullpostY,2,mean)
# Put in matrix consistent with gjam output for comparisons
yMeanMat <- matrix(yMean_gjam, nrow=dim(ObsY)[1], ncol=dim(ObsY)[2],byrow = F)
head(yMeanMat)
## Plot observed vs. predicted (marginal predictions)
# Predicted values of Y
# Posterior means
yMu  <- jdm1$prediction$ypredMu 
yMuv <- as.vector(yMu)
# Observed values of Y
ObsY <- jdm1$inputs$y
ObsYv <- as.vector(ObsY)


plot(ObsYv, yMuv)
# Posterior medians
points(ObsYv, yMedian_gjam, col='red')
# points(ObsYv, yMean_gjam, col='green')
abline(0,1)

# plot(yMu[,1], ObsY[,1])
# abline(0,1)
# 
# plot(yMu[,3], ObsY[,3])
# abline(0,1)

######################################
######################################
# MARGINAL prediction can also be done usig gjamPredict
# Obtain full posterior predictive distributions using FULL=T
xdata <- xdat
newdata   <- list(xdata = xdata, nsim = 1000 )
p1 <- gjamPredict(jdm1, newdata = newdata, FULL=T)
# str(p1)

plot(ObsY, p1$sdList$yMu)
abline(0,1)

### p1$ychains has the posterior predictive distributions for each observation and species

# Grab posterior samples for a species
# spp1 <- p1$ychains[,1:1343]
# spp1Mean <- apply(spp1,2,mean)
# hist(spp1Mean)
# spp1Median <- apply(spp1,2,median)
# hist(spp1Median)

# Obtain posterior medians for all obs
yMedian <- apply(p1$ychains,2,median)
length(yMedian)

ObsYv <- as.vector(ObsY)

plot(ObsYv, yMedian)
abline(0,1)


# Change vector of predictions to matrix for plotting by species
yMedianMat <- matrix(yMedian, nrow=dim(ObsY)[1], ncol=dim(ObsY)[2],byrow = F)

# black bullhead
plot(yMu[,1], ObsY[,1])
points(yMedianMat[,1],ObsY[,1], col='red')
abline(0,1)

# Walleye
plot(yMu[,12], ObsY[,12])
points(yMedianMat[,12],ObsY[,12], col='red')
abline(0,1)

# Yellow perch
plot(yMu[,15], ObsY[,15])
points(yMedianMat[,15],ObsY[,15], col='red')
abline(0,1)




# prediction$presence gives you the probability of presence for the species/location. 
# Therefore, 1-prediction$presence gives the probability of a 0
probZero <- 1-jdm1$prediction$presence
head(probZero)
# Replace post means with zero if prob zero greater than threshold
yMu[probZero >0.7] <- 0
head(yMu)
apply(probZero,2,range)
plot(yMu, ObsY)
abline(0,1)

# # Walleye
# plot(yMu[,12], ObsY[,12])
# abline(0,1)
# 
# # Perch
# plot(yMu[,14], ObsY[,14])
# abline(0,1)
# 
# plot(yMu[,3], ObsY[,3])
# abline(0,1)

# plot(probZero[,3],ObsY[,3])


###### Plot average CPE  after removing the 0s to see fits. 
# Ave CPE is computed as Ymu divided by the probability of presence. 
# Predicted values of Y
yMu  <- jdm1$prediction$ypredMu 
head(yMu)
summary(yMu)
# Remove zero preds
yMu[yMu ==0] <- NA

yMuAve <- yMu/jdm1$prediction$presence
head(yMuAve)

plot(yMuAve, ObsY)
abline(0,1)

plot(yMuAve[,3], ObsY[,3])
abline(0,1)



# Posterior means of beta
jdm1$parameters$betaMu

# Posterior means of sigma
jdm1$parameters$sigMu
dim(jdm1$parameters$sigMu)

# Posterior means of Rho
jdm1$parameters$corMu

head(jdm1$chains$bgibbs) #these are all the betas
BetaOut <- jdm1$chains$bgibbs
dim(BetaOut)

# apply(BetaOut,2,mean)

head(jdm1$chains$sgibbs) #these are all the elements of Sigma
SigmaOut <-jdm1$chains$sgibbs
dim(SigmaOut)

# # Make upper triangular matrix set to zero 
# RhoCILLagos[upper.tri(RhoCILLagos)] <- NA
# # remove diagonal
# RhoCILLagos[upper.tri(RhoCILLagos,diag=TRUE)] <- NA


saveRDS(BetaOut, file="BetaOut.rds")
saveRDS(SigmaOut, file="SigmaOut.rds")
saveRDS(jdm1, file="gjamOUT1.rds")



