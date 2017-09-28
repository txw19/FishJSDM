# rm(list=ls())

# Read in saved data and compute stats
##### Residual correlations
Rho <- readRDS(file="ResRho.rds")

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


########################################################
########################### Environmental Correlations

ERho <- readRDS(file="EnvRho.rds")

ERhoMean <- apply(ERho,c(3,2),mean)
# Make upper triangular matrix set to zero 
ERhoMean[upper.tri(ERhoMean)] <- NA
# remove diagonal
ERhoMean[upper.tri(ERhoMean,diag=TRUE)] <- NA
ERMean <- as.numeric(t(ERhoMean))
ERMean <- ERMean[!is.na(ERMean)]

# Lower CI
ERhoCIL <- apply(ERho,c(3,2),quantile, 0.025)
# Make upper triangular matrix set to zero 
ERhoCIL[upper.tri(ERhoCIL)] <- NA
# remove diagonal
ERhoCIL[upper.tri(ERhoCIL,diag=TRUE)] <- NA
ERCIL <- as.numeric(t(ERhoCIL))
ERCIL <- ERCIL[!is.na(ERCIL)]

# Upper CI
ERhoCIU <- apply(ERho,c(3,2),quantile, 0.975)
# Make upper triangular matrix set to zero 
ERhoCIU[upper.tri(ERhoCIU)] <- NA
# remove diagonal
ERhoCIU[upper.tri(ERhoCIU,diag=TRUE)] <- NA
ERCIU <- as.numeric(t(ERhoCIU))
ERCIU <- ERCIU[!is.na(ERCIU)]



# sppnames <- c("black bullhead","black crappie","bowfin","brown bullhead",
#               "common carp","golden shiner","largemouth bass","northern pike",
#               "rock bass","smallmouth bass","cisco","walleye",
#               "white sucker","yellow bullhead","yellow perch","sunfish")

# Finding example correlations
min(RMean)
ERhoMean
which(ERMean==min(ERMean))
RMean[118]

which(RMean==max(RMean))
ERMean[112]

RMean[3]
ERMean[3]

#####################################################
########### PLOT ####################################
#####################################################

res <- 6
name_figure <- "CorrPanel.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...


nf <- layout(matrix(c(1:1), ncol=1, nrow = 1,byrow=TRUE), respect=T)
# , respect=T
layout.show(nf)
# par(mar = c(1, 0.5, 0, 0) + 0.1,oma=c(1.5,1.0,0,0))
par(oma=c(2.1,2.0,0,0),mai=c(0.0,0.00,0.0,0))


def.par <- par(no.readonly = TRUE)     

size.labels = 1
size.text = 0.9

x.label <- 'Environmental correlation'
y.label <- 'Residual correlation'

# Plot
plot(ERMean, RMean,xlim=c(-1,1),ylim=c(-1,1),  
     axes=F, xlab='',ylab='',type='n')
  
axis(side=2,at=c(-1, -0.5, 0, 0.5, 1),labels=c(-1, -0.5, 0, 0.5, 1),tck= -0.005, mgp=c(0,0.3,0), las=1)

axis(side=1,cex.axis=size.text , mgp=c(1,0,0),tck= -0.005, labels=T,at=c(-1, -0.5, 0, 0.5, 1) ) 

abline(v=0, col="black")
abline(h=0, col="black")

# 95% CIs for censored analysis
segments(x0=ERCIU, x1=ERCIL,
         y0=RMean, y1=RMean, col="gray",lwd=1)

segments(x0=ERMean, x1=ERMean,
         y0=RCIU, y1=RCIL, col="gray",lwd=1)

points(ERMean, RMean, col="black" ,cex=0.8, pch=16)

# Add x- and y-axis lables
mtext(x.label, line = 1, side = 1, cex = 1, outer=T, adj=0.5)
mtext(y.label, line = 1, side = 2, cex = 1, outer=T,adj=0.5)
box()

# arrows(0.8,0.6,0.92,0.47, length=0.15)
# text(0.8, 0.62, "Largemouth and sunfish", cex=0.7)
# 
# arrows(-0.94,-0.45,-0.94,-0.20, length=0.15)
# text(-0.8, -0.5, "White sucker and sunfish", cex=0.7)

par(def.par)
dev.off()

