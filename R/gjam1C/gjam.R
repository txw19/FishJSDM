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

fishP=fread("catch_and_effort_most_recent_with_predictors.csv")
# fishP=fread("cpue_JDM.csv")
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


# Remove column
dat[,dowlknum:=NULL]

# Remove lakes with missing Secchi
dat <- dat[!is.na(dat$Secchi.lake.mean),]

# Find some "outliers"
which(dat$black.crappie>2600)
which(dat$brown.bullhead>500)
which(dat$golden.shiner>300)
# Remove outliers
dat <- dat[dat$black.crappie < 2600 & dat$brown.bullhead < 500 & dat$golden.shiner < 300]
dim(dat)

###########################
# Read in lake lat/longs
lls <- fread("mn_lake_list.csv")
head(lls)
dim(lls)
# Grab lat/long columns
lls2 <- lls[, .(DOW, LAKE_CENTER_LAT_DD5,LAKE_CENTER_LONG_DD5)]
head(lls2)
dim(lls2)

setkey(dat, DOW)
setkey(lls2, DOW)
dat <- merge(dat, lls2)
dat[,.N]
head(dat)

# Erin data-add predictions
Edat <- cbind(dat,yMu)
dim(Edat)

head(Edat)
# write.csv(Edat, "FinaldatForE.csv", row.names = F)

# write.csv(dat, "Finaldat.csv", row.names = F)



summary(dat)


xdat <- dat[,23:27]
ydat <- dat[,4:19]

xdat <- as.data.frame(xdat)
ydat <- as.data.frame(ydat)

# Effort
ef <- as.numeric(dat$EFFORT)

ef  <- list( columns = 1:16, values = dat$EFFORT)


# Some gjam options for modelList:
# holdoutN = 0, number of observations to hold out for out-of-sample prediction.
# holdoutIndex =  numeric(0), numeric vector of observations (row numbers) to holdout for out-of-sample prediction
# ,holdoutN=200
start.time = Sys.time()  # Start timer: 56 minute run time

ml <-list(PREDICTX = F, ng=70000,burnin=60000,typeNames=rep("DA",16), effort=ef) 
ml$FULL=T

# jdm1 = gjam(~ log_area + log_depth + Secchi.lake.mean + mean.gdd + alkalinity, 
#             xdata=xdat, ydata=ydat,
#             modelList=ml)
# I(temp^2)

jdm1 = gjam(~ log_area + log_depth + Secchi.lake.mean + mean.gdd + alkalinity + 
              I(log_area^2) + I(log_depth^2) + I(Secchi.lake.mean^2)+ I(mean.gdd^2) + I(alkalinity^2),
          xdata=xdat, ydata=ydat,
          modelList=ml)

end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='')

# jdm1 <- readRDS(file="gjamOUT1.rds")
str(jdm1)
summary(jdm1)
names(jdm1)

# jdm1$parameters$corMu

# Plot gjam output
plotPars <- list(GRIDPLOTS=T,PLOTTALLy=F,SAVEPLOTS=T,SMALLPLOTS = F) 
fit <- gjamPlot(jdm1, plotPars)

# in-sample posterior predictions
fullpostY <- jdm1$chains$ygibbs
# dim(fullpostY)
# # Calculate posterior medians
yMedian_gjam <- apply(fullpostY,2,median)
# length(yMedian_gjam)

# Posterior means
yMu  <- jdm1$prediction$ypredMu 
yMuv <- as.vector(yMu)
# Observed values of Y
ObsY <- jdm1$inputs$y
ObsYv <- as.vector(ObsY)

# Plot observed vs. posterior mean
plot(ObsYv, yMuv)
# Overlay posterior medians
# points(ObsYv, yMedian_gjam, col='red')
# points(ObsYv, yMean_gjam, col='green')
abline(0,1)

# In-sample predictions, i.e., "xdata" is not provided in gjamPredict
gjamPredict(jdm1, y2plot = colnames(ydat))
# Northern Pike
gjamPredict(jdm1, y2plot = colnames(ydat)[8])
# Walleye
gjamPredict(jdm1, y2plot = colnames(ydat)[12])
gjamPredict(jdm1, y2plot = colnames(ydat)[13])

# Posterior predictions for each species
sppPred <- array(fullpostY, dim=c(dim(jdm1$chains$bgibbs)[1],dim(xdat)[1],dim(jdm1$prediction$presence)[2]))
hist(sppPred[,,1])


# Figure of predictions for each species
# sppnames <- c("black bullhead","black crappie","bowfin","brown bullhead",
#               "common carp","golden shiner","largemouth bass","northern pike",
#               "rock bass","smallmouth bass","cisco","walleye",
#               "white sucker","yellow bullhead","yellow perch","sunfish")
res <- 6
name_figure <- "ObsPredBySPP_Panel.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...

nf <- layout(matrix(c(1:16), ncol=4, nrow = 4,byrow=TRUE), respect=F)
layout.show(nf)

par(mar = c(1, 0, 0, 0) + 0.1,oma=c(2,2.5,0,0.5))

for(i in 1:16){
  s1 <- apply(sppPred[,,i],2,median)
  s2 <- apply(sppPred[,,i],2,mean)
  plot(ObsY[,i],s2,axes=F, xlab='',ylab='',type='p', cex=0.5)
  # points(ObsY[,i],s2, col="red", cex=0.5)
  axis(side=1,cex.axis=0.8, mgp=c(1,0,0),tck= -0.01)
  axis(side=2,cex.axis=0.8, tck= -0.005, mgp=c(0,0.3,0), las=1)
  box()
  abline(0,1)
  mtext("Obs", line = 0.8, side = 1, cex = 1, outer=T, adj=0.5)
  mtext("pred", line = 1, side = 2, cex = 1, outer=T, adj=0.5)
}
par(def.par)
dev.off()



######################################
######################################
# Out-of-sample prediction can also be done usig gjamPredict, specifying xdata
# Obtain full posterior predictive distributions using FULL=T
xdata <- xdat
newdata   <- list(xdata = xdata, nsim = 1000 )
p1 <- gjamPredict(jdm1, newdata = newdata, FULL=T)
# str(p1)

plot(ObsY, p1$sdList$yMu)
abline(0,1)

### p1$ychains has the posterior predictive distributions for each observation and species



# head(jdm1$chains$bgibbs) #these are all the betas
BetaOut <- jdm1$chains$bgibbs
dim(BetaOut)
hist(BetaOut[,3])
# apply(BetaOut,2,mean)

# head(jdm1$chains$sgibbs) #these are all the elements of Sigma
SigmaOut <-jdm1$chains$sgibbs
dim(SigmaOut)


saveRDS(BetaOut, file="BetaOut.rds")
saveRDS(SigmaOut, file="SigmaOut.rds")
# saveRDS(jdm1, file="gjamOUT1.rds")



