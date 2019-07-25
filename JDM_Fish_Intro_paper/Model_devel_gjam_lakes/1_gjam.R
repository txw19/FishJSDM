# rm(list=ls())
# load the package
library(devtools)


# library(Hmsc)
library(gjam)
library(dplyr)
library(coda)
library(xtable)
library(ggplot2)
library(data.table)
library(ggsn) # add N arrow
library(rgdal)
library(maps)
library(mapdata)
library(ggmap)
library(sf)
library(ggrepel) # for labeling map
library(doBy)
library(car) # for the logit function
library(sp)
library(lubridate)
library(colorspace)
library(corrplot)
library(circlize)
library(psych)
library(formattable)
library(knitr)
library(tidyr)
library(Hmisc)

# read in data
fishP=fread("lake_predictors_fish_PA_for_HMSC.csv")

dim(fishP) 
fishP
summary(fishP)

# Remove lakes with missing (NA) predictors
fishP <- fishP[!is.na(littoral.zone)]
fishP <- fishP[!is.na(mean.gdd)]
fishP <- fishP[!is.na(proportion_disturbed_ws)]
summary(fishP)
dim(fishP)

# Export final data
write.csv(fishP, "finalDat.csv", row.names = F)
# # #Remove species with low detections (present at <10 sites)
# # fishP[,c(9,14,15,16,17,19,21,23,26,29,30,31) ]<- NULL
# fishP[,c(9,14,11,15,16,17,18,19,20,21,23,25,26,29,30,31) ]<- NULL

# # Read in and merge covariates
# covs <- fread('lake_predictors_fish_PA_for_HMSC.csv')
# dim(covs)
# 
# fishP <-plyr::join(fishP, covs, by=c("DOW"), type='left', match='all')
# dim(fishP)
# head(fishP)
# summary(fishP)
##########################################
### Get data in order for hmsc
# Species matrix, a site * species matrix
ydat <- fishP[,8:24]
ydat <- as.data.frame(ydat)



# #Determining co-linearity 
# MyVar<-fishP[,c(21:27)]
# pairs.panels(MyVar, method = "pearson", hist.col = "#00AFBB", density = TRUE, ellipses = FALSE)

# Covariate matrix
xdat<- fishP[,2:7]
xdat <- data.table(xdat)
head(xdat)
xdat$area.hectares <- log(xdat$area.hectares)
summary(xdat)
# # Change all covariate data columns to numeric
# changeCols <- colnames(xdat)[2:6]
# xdat[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols]

xdat <- as.data.frame(xdat)
str(xdat)
head(xdat)
# pairs.panels(xmat, method = "pearson", hist.col = "#00AFBB", density = TRUE, ellipses = FALSE)


# REDUCT = F in modelList overrides automatic dimension reduction.
ml <-list(FULL=T, PREDICTX = F, ng=90000,burnin=70000,typeNames=rep("PA",dim(ydat)[2]), REDUCT=F) 


jdm1 = gjam(~ area.hectares + max_depth_m + mean.gdd + proportion_disturbed_ws + littoral.zone + secchi.m,
            xdata=xdat, ydata=ydat,
            modelList=ml)

saveRDS(jdm1, file="gjamOUT1.rds")
# jdm1 <- readRDS(file="gjamOUT1.rds")
str(jdm1)
summary(jdm1)
names(jdm1)

# Grab residual correlation
# jdm1$parameters$corMu
# 
# jdm1$parameter$rndEff
# jdm1$parameter$sigMu

# Plot gjam output
plotPars <- list(GRIDPLOTS=T,PLOTTALLy=F,SAVEPLOTS=T,SMALLPLOTS = F) 
fit <- gjamPlot(jdm1, plotPars)

#### Species richness (responses predicted > 0)
rich <- jdm1$prediction$richness
head(rich)
dim(rich)


########## marginal prediction
newdata <- list(xdata = xdat, nsim=10000)
tmp     <- gjamPredict(jdm1, newdata=newdata)
full    <- tmp$sdList$yMu

# Calculate marginal AUC for each species (larger is better)
AUC <- numeric()
for(i in 1:dim(ydat)[2]){
  AUC[i] <- somers2(full[,i], ydat[,i])[1]
}
AUC


## Conditional predictions
############## For all species
# Grab species names for looping through species 
cnames <- colnames(ydat) 
cond.preds <- matrix(NA, nrow = 2, ncol=length(cnames))
dim(cond.preds)

#right now output has coldwaterspp, will change maybe on next run
cnames[17]="coldwaterspp"
colnames(ydat)[17]="coldwaterspp"


for(j in 1:length(cnames)){ 
  yc <- ydat[ , !(colnames(ydat) %in% cnames[j])] # Loop through each species for conditional predictions
  #set species of interest to 0, 
 
  new <- list(ydataCond = yc, xdata=xdat, nsim=10000)   # cond on obs P/A data 
  preds <- gjamPredict(output = jdm1, newdata = new) 
  condy  <- preds$sdList$yMu
  cond.preds[2,j] <- somers2(condy[,(colnames(condy) %in% cnames[j])], ydat[,(colnames(condy) %in% cnames[j])])[1]
  cond.preds[1,j] <- somers2(full[,(colnames(condy) %in% cnames[j])], ydat[,(colnames(condy) %in% cnames[j])])[1]
  print(j)
}

colnames(cond.preds) <- cnames
rownames(cond.preds) <- c("Marginal AUC", "Conditional AUC")

write.csv(cond.preds, "Marginal_conditional_AUCs.csv")

# MCMC samples of Sigma
SigmaOut <-jdm1$chains$sgibbs
str(SigmaOut)
saveRDS(SigmaOut, file="SigmaOut.rds")

## Convert Sigma to Corr
# SigmaOut.rds contains the lower triangle of the variance covariance matrix, e.g., 136 columns with 16 species
Sigma <- readRDS(file="SigmaOut.rds")
str(Sigma)
dim(Sigma)

postSigMeans <- apply(Sigma,2,mean)
length(postSigMeans)

S <- diag(dim(ydat)[2])
S[lower.tri(S, diag=TRUE)] <- postSigMeans

# Convert S (half matrix) to full matrix
X <- diag(dim(ydat)[2])
X[lower.tri(X, diag=TRUE)] <- postSigMeans
X <- X + t(X) - diag(diag(X))


# Convert Sigma to Rho
Rho <- array(NA, dim=c(dim(Sigma)[1], dim(X)[1], dim(X)[1]))
dim(Rho)
for(samp in 1:dim(Sigma)[1]){
  # Grab mcmc sample from Sigma
  sigsamp <- Sigma[samp,]
  # Convert S (half matrix) to full matrix
  Stemp <- diag(dim(ydat)[2]) 
  Stemp[lower.tri(Stemp, diag=TRUE)] <- sigsamp 
  Stemp <- Stemp + t(Stemp) - diag(diag(Stemp)) 
  # Calculate correlation matrix from Sigma
  Rho[samp, , ] <- cov2cor(Stemp)
}


dim(Rho)
# Posterior mean
RhoMean <- apply(Rho,c(3,2),mean)
write.csv(RhoMean,'resid.correlation.matrix.mn.lakes.csv')

# Lower CI
RhoCIL <- apply(Rho,c(3,2),quantile, 0.025)

# Upper CI
RhoCIU <- apply(Rho,c(3,2),quantile, 0.975)

# Use this to plot only significant cors
sigCors <- RhoCIL * RhoCIU > 0


########### Residual correlations
species.names<-colnames(ydat)
# Posterior means of Rho
# jdm1$parameters$corMu
# colnames(ydat)
colnames(RhoMean) <- species.names
rownames(RhoMean) <- species.names

colMat <- matrix(NA, nrow = nrow(RhoMean), ncol = ncol(RhoMean))
colMat[which(RhoMean > 0 & sigCors =='TRUE', arr.ind = TRUE)] <- "blue"
colMat[which(RhoMean < 0 & sigCors =='TRUE', arr.ind = TRUE)] <- "red"

res<-6
tiff(filename = 'resid_corr_mn_lakes2.tiff', height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE) 
chordDiagram(RhoMean, symmetric = FALSE,
             annotationTrack = c("grid"), grid.col = "grey",
             col = colMat,
             annotationTrackHeight = c(0.01, 0.01),
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(RhoMean))))))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

par(def.par)
dev.off()

# Posterior means
# sigma
# SigmaOut <-jdm1$chains$sgibbs
BetaOut <- jdm1$chains$bgibbs
saveRDS(BetaOut, file="BetaOut.rds")
# saveRDS(SigmaOut, file=paste("SigmaOut",i,".rds",sep=''))
# saveRDS(jdm1, file=paste("gjamOut",i,".rds",sep=''))


