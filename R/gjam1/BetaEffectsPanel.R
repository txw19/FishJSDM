# rm(list=ls())


# Read in saved data and compute stats
betaout <- readRDS(file="BetaOut.rds")

# Remove intercepts from Beta

int.names <- c("black.bullhead_intercept","black.crappie_intercept","bowfin.(dogfish)_intercept","brown.bullhead_intercept",
               "common.carp_intercept","golden.shiner_intercept","largemouth.bass_intercept","northern.pike_intercept",
               "rock.bass_intercept","smallmouth.bass_intercept","tullibee.(cisco)_intercept","walleye_intercept",
               "white.sucker_intercept","yellow.bullhead_intercept","yellow.perch_intercept","sunfish_intercept")

betaout <- betaout[, !colnames(betaout) %in% int.names]
dim(betaout)

betaMeans <- apply(betaout,2,mean)
betaCRI <- apply(betaout,2,quantile, c(0.025, 0.975))


betaMeans <- t(betaMeans)
dim(betaMeans)


sppnames <- c("black bullhead","black crappie","bowfin","brown bullhead",
              "common carp","golden shiner","largemouth bass","northern pike",
              "rock bass","smallmouth bass","cisco","walleye",
              "white sucker","yellow bullhead","yellow perch","sunfish")

#####################################################
########### PLOT ####################################
#####################################################

covariates <- c("Area", "Depth", "Secchi","GDD", "Alkalinity")

res <- 6
name_figure <- "FishEffectsPanel.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...

nf <- layout(matrix(c(1:16), ncol=8, nrow = 2,byrow=TRUE), respect=F)
# , respect=T
layout.show(nf)
# par(mar = c(1, 0.5, 0, 0) + 0.1,oma=c(1.5,1.0,0,0))
# par(mar=c(0.0,0.0,0.0,0.0),oma=c(1,2.5,0,0),mai=c(0.0,0.03,0.03,0) )
par(mar = c(1, 0, 0, 0) + 0.1,oma=c(2,2.5,0,0.5))

size.labels = 1
size.text = 0.9

x.label <- 'Estimated effect'
y.label <- expression(paste(beta ,' parameters'))

# Posterior means and CIs for all parameters
Plot.data <- rbind(betaMeans,betaCRI)

# 16 = number of species
rows <- 1:(dim(Plot.data)[2]/16)

plotting.region <- range(betaCRI)

### axis label options
spc <- 0.55
lab <- 1:63
cex <- 0.7
adj <- 0
###

# Dimension of Plot.data, need to plot e.g., 1:3 for bluegill, 4:6, 7:9, 10:12, 13:15, 16:18 for 
# each species and 3 covariates (6 * 3 = 18)
# Thus, c(1,4,7,10,13,16) in for loop corresponds to first column to plot for each species.
# We then add 4 to i to get the last column, e.g., Plot.data[1,i:(i+2)]


for(i in c(1,6,11,16,21,26,31,36,41,46,51,56,61,66,71,76)){
# Color for plotting
Plot.color <- as.numeric(betaCRI[1,i:(i+4)] * betaCRI[2,i:(i+4)] > 0 )
colorP <- rep("blue", length(rows))
colorP[Plot.color > 0] <- "black"

plotting.region <- range(betaCRI[,i:(i+4)])

plot(c(plotting.region[1]-adj, plotting.region[2]), c(0.5,length(rows)+0.5), 
     axes=F, xlab='',ylab='',type='n')
  

if(i==1|i==41){axis(side=2,at=c(1:length(rows)),labels=covariates,tck= 0.001, mgp=c(0,0.3,0)) }
axis(side=1,cex.axis=size.text , mgp=c(1,0,0),tck= -0.01)


# 95% CIs for censored analysis
segments(x0=Plot.data[2,i:(i+4)], x1=Plot.data[3,i:(i+4)],
         y0=1:length(rows), y1=1:length(rows), col=colorP,lwd=1, lty=1)

points(Plot.data[1,i:(i+4)], 1:length(rows), col=colorP ,cex=1, pch=16)

abline(v=0, col="gray")

mtext(x.label, line = 0.8, side = 1, cex = 1, outer=T, adj=0.5)
mtext(y.label, line = 1, side = 2, cex = 1, outer=T)

# for(i in c(1,6,11,16,21,26,31,36,41,46,51,56,61,66,71,76)){

if(i==1){ text(plotting.region[1]+1.4, 5.5, sppnames[1], cex=0.8) }
if(i==6){ text(plotting.region[1]+0.7, 5.5, sppnames[2], cex=0.8) }
if(i==11){ text(plotting.region[1]+0.2, 5.5, sppnames[3], cex=0.8) }
if(i==16){ text(plotting.region[1]+0.3, 5.5, sppnames[4], cex=0.8) }
if(i==21){ text(plotting.region[1]+0.9, 5.5, sppnames[5], cex=0.8) }
if(i==26){ text(plotting.region[1]+0.4, 5.5, sppnames[6], cex=0.8) }
if(i==31){ text(plotting.region[1]+0.2, 5.5, sppnames[7], cex=0.8) }
if(i==36){ text(plotting.region[1]+0.15, 5.5, sppnames[8], cex=0.8) }

if(i==41){ text(plotting.region[1]+0.6, 5.5, sppnames[9], cex=0.8) }
if(i==46){ text(plotting.region[1]+0.4, 5.5, sppnames[10], cex=0.8) }
if(i==51){ text(plotting.region[1]+1.5, 5.5, sppnames[11], cex=0.8) }
if(i==56){ text(plotting.region[1]+0.5, 5.5, sppnames[12], cex=0.8) }
if(i==61){ text(plotting.region[1]+0.5, 5.5, sppnames[13], cex=0.8) }
if(i==66){ text(plotting.region[1]+0.4, 5.5, sppnames[14], cex=0.8) }
if(i==71){ text(plotting.region[1]+0.8, 5.5, sppnames[15], cex=0.8) }
if(i==76){ text(plotting.region[1]+0.4, 5.5, sppnames[16], cex=0.8) }

# abline(h=0)
box()
} # end loop
par(def.par)
dev.off()

