# rm(list=ls())
library(gjam)
library(data.table)

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

xdat <- dat[,20:24]
ydat <- dat[,2:17]

# Some gjam options for modelList:
# holdoutN = 0, number of observations to hold out for out-of-sample prediction.
# holdoutIndex =  numeric(0), numeric vector of observations (row numbers) to holdout for out-of-sample prediction

# -typeNames = "CA" means continuous abundance (meaning greater than 0).
start.time = Sys.time()  # Start timer 

jdm1 = gjam(~ log_area + log_depth + Secchi.lake.mean + mean.gdd + alkalinity, xdata=xdat, ydata=ydat,
          modelList=list(PREDICTX = F, ng=90000,burnin=50000,typeNames=rep("CA",16)))

end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

summary(jdm1)

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



