# rm(list=ls())


# Read in saved data and compute stats
betaout <- readRDS(file="BetaOut.rds")
betaMeans <- apply(betaout,2,mean)
betaCRI <- apply(betaout,2,quantile, c(0.025, 0.975))


betaMeans <- t(betaMeans)
dim(betaMeans)

# remove intercept estimates for plotting
betaMeans <- betaMeans[-c(1,5,9,13,17,21)]
betaCRI <- betaCRI[, -c(1,5,9,13,17,21)]

#####################################################
########### PLOT ####################################
#####################################################

covariates <- c("Depth", "GDD", "Secchi")

res <- 6
name_figure <- "FishEffectsPanel.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...

nf <- layout(matrix(c(1:6), ncol=6, nrow = 1,byrow=TRUE), respect=F)
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

# 6 = number of species
rows <- 1:(dim(Plot.data)[2]/6)

plotting.region <- range(betaCRI)

### axis label options
spc <- 0.55
lab <- 1:63
cex <- 0.7
adj <- 0
###

species <- c("Bluegill", "Northern Pike", "Walleye", "White Sucker", "Yellow perch", "Largemouth bass" )

# Dimension of Plot.data, need to plot e.g., 1:3 for bluegill, 4:6, 7:9, 10:12, 13:15, 16:18 for 
# each species and 3 covariates (6 * 3 = 18)
# Thus, c(1,4,7,10,13,16) in for loop corresponds to first column to plot for each species.
# We then add 2 to i to get the last column, e.g., Plot.data[1,i:(i+2)]

for(i in c(1,4,7,10,13,16)){
# Color for plotting
Plot.color <- as.numeric(betaCRI[1,i:(i+2)] * betaCRI[2,i:(i+2)] > 0 )
colorP <- rep("blue", length(rows))
colorP[Plot.color > 0] <- "black"

plotting.region <- range(betaCRI[,i:(i+2)])

plot(c(plotting.region[1]-adj, plotting.region[2]), c(0.5,length(rows)+0.5), 
     axes=F, xlab='',ylab='',type='n')
  

if(i==1){axis(side=2,at=c(1:length(rows)),labels=covariates,tck= 0.001, mgp=c(0,0.3,0)) }
axis(side=1,cex.axis=size.text , mgp=c(1,0,0),tck= -0.01)

# if( i <=10){
#   axis(side=1,cex.axis=size.text , mgp=c(1,0,0),tck= -0.01, labels=F ) 
# } else {
#   axis(side=1,cex.axis=size.text , mgp=c(1,0,0),tck= -0.01)
# }	
# if(i==1){ text(-8,1:length(rows),srt = 0, adj =adj,labels = covariates, xpd = TRUE,cex=cex) }


# 95% CIs for censored analysis
segments(x0=Plot.data[2,i:(i+2)], x1=Plot.data[3,i:(i+2)],
         y0=1:length(rows), y1=1:length(rows), col=colorP,lwd=1, lty=1)

points(Plot.data[1,i:(i+2)], 1:length(rows), col=colorP ,cex=1, pch=16)

abline(v=0, col="gray")

mtext(x.label, line = 0.8, side = 1, cex = 1, outer=T, adj=0.5)
mtext(y.label, line = 0.9, side = 2, cex = 1, outer=T)

if(i==1){ text(plotting.region[1]+5, 3.5, species[1]) }
if(i==4){ text(plotting.region[1]+0.7, 3.5, species[2]) }
if(i==7){ text(plotting.region[1]+0.7, 3.5, species[3]) }
if(i==10){ text(plotting.region[1]+0.9, 3.5, species[4]) }
if(i==13){ text(plotting.region[1]+4.5, 3.5, species[5]) }
if(i==16){ text(plotting.region[1]+0.37, 3.5, species[6]) }


# abline(h=0)
box()
} # end loop
par(def.par)
dev.off()

