# rm(list=ls())
library(ggplot2)
library(reshape2)
library(ggrepel)
library(dplyr)

setwd("~/MN fish trends/Joint Distribution Model/JDM_Fish_Intro_paper/All_species/")
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

finalDat=read.csv("finalDat.csv")

 sppnames <- colnames(finalDat)[6:21]


colnames(ERhoMean)=sppnames
rownames(ERhoMean)=sppnames 


colnames(ERhoCIU)=sppnames
rownames(ERhoCIU)=sppnames 

colnames(ERhoCIL)=sppnames
rownames(ERhoCIL)=sppnames 

corr.data.CIU=melt(ERhoCIU, value.name = "Env.cor.CIU", varnames = c("Species1", "Species2"))
corr.data.CIL=melt(ERhoCIL, value.name = "Env.cor.CIL", varnames = c("Species1", "Species2"))
corr.data.mean=melt(ERhoMean, value.name = "Env.cor", varnames = c("Species1", "Species2"))
 
corr.data=merge(corr.data.mean, corr.data.CIL, by=c("Species1", "Species2"))
corr.data=merge(corr.data, corr.data.CIU, by=c("Species1", "Species2"))

##################################

colnames(RhoMean)=sppnames
rownames(RhoMean)=sppnames 

colnames(RhoCIU)=sppnames
rownames(RhoCIU)=sppnames 

colnames(RhoCIL)=sppnames
rownames(RhoCIL)=sppnames 

corr.data2.CIU=melt(RhoCIU, value.name = "Res.cor.CIU", varnames = c("Species1", "Species2"))
corr.data2.CIL=melt(RhoCIL, value.name = "Res.cor.CIL", varnames = c("Species1", "Species2"))
corr.data2.mean=melt(RhoMean, value.name = "Res.cor", varnames = c("Species1", "Species2"))
 
corr.data2=merge(corr.data2.mean, corr.data2.CIL, by=c("Species1", "Species2"))
corr.data2=merge(corr.data2, corr.data2.CIU, by=c("Species1", "Species2"))


corr.data=merge(corr.data, corr.data2, by=c("Species1", "Species2"))
corr.data=corr.data[complete.cases(corr.data),]


#if going to facet, need to repeat so all species dots appear in all species panels, but switching species 1 and species 2


corr.data.doubled=corr.data
corr.data.doubled$Species2=corr.data$Species1
corr.data.doubled$Species1=corr.data$Species2


corr.data.plot=rbind(corr.data.doubled, corr.data)

windows()
ggplot(corr.data.plot, aes(Env.cor, Res.cor, colour=Species1))+geom_hline(yintercept = 0, alpha=0.5, colour="grey")+geom_vline(xintercept=0, alpha=0.5, colour="grey")+geom_point()+geom_text_repel(aes(label=Species1), size=3) +facet_wrap(~Species2)+theme_bw()+xlab("Shared environmental correlation")+ylab("Residual correlation")+theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank())
ggsave("environmental_and_residual_correlations_by_species_PA_streams.png", height=12, width=10, units="in")

#just CIs that don't overlap zero
corr.data.plot$plot.env=0
corr.data.plot$plot.env[sign(corr.data.plot$Env.cor.CIL)==sign(corr.data.plot$Env.cor.CIU)]=1

corr.data.plot$plot.res=0
corr.data.plot$plot.res[sign(corr.data.plot$Res.cor.CIL)==sign(corr.data.plot$Res.cor.CIU)]=1

corr.data.plot2=subset(corr.data.plot, plot.res==1 |plot.env==1)

windows()
ggplot(corr.data.plot2, aes(Env.cor, Res.cor, colour=Species1))+geom_hline(yintercept = 0, alpha=0.5, colour="grey")+geom_vline(xintercept=0, alpha=0.5, colour="grey")+geom_point()+geom_text_repel(aes(label=Species1), size=3) +facet_wrap(~Species2)+theme_bw()+xlab("Shared environmental correlation")+ylab("Residual correlation")+theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank())+geom_errorbarh(aes(xmin=Env.cor.CIL, xmax=Env.cor.CIU), alpha=.25)+geom_errorbar(aes(ymin=Res.cor.CIL, ymax=Res.cor.CIU), alpha=.25)
ggsave("environmental_and_residual_correlations_by_species_PA_significant_streams.png", height=12, width=10, units="in")


#can we make an environmental chord plot?
########### Env correlations

colnames(ERhoMean) <- sppnames
rownames(ERhoMean) <- sppnames

colMat <- matrix(NA, nrow = nrow(ERhoMean), ncol = ncol(ERhoMean))
# Lower CI
ERhoCIL <- apply(ERho,c(3,2),quantile, 0.025)

# Upper CI
ERhoCIU <- apply(ERho,c(3,2),quantile, 0.975)

# Use this to plot only significant cors
sigCors <- ERhoCIL * ERhoCIU > 0

colMat[which(ERhoMean > 0 & sigCors =='TRUE', arr.ind = TRUE)] <- "blue"
colMat[which(ERhoMean < 0 & sigCors =='TRUE', arr.ind = TRUE)] <- "red"

res<-6
tiff(filename = 'env_corr_streams.tiff', height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE) 
chordDiagram(ERhoMean, symmetric = FALSE,
             annotationTrack = c("grid"), grid.col = "grey",
             col = colMat,
             annotationTrackHeight = c(0.01, 0.01),
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(ERhoMean))))))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

par(def.par)
dev.off()


#plot AUC
AUC=read.csv("Marginal_conditional_AUCs.csv")
head(AUC)
colnames(AUC)[1]="type"

AUC.long=melt(AUC, variable.name = "species")

# The arguments to spread():
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
AUC.wide <- spread(AUC.long, type, value)
colnames(AUC.wide)=c("species", "conditional.auc", "marginal.auc")

#arbitrary thresholds for accurate vs inaccurate models
thresholds=data.frame("accurate"=0.9, "inaccurate"=0.7, edge1=Inf, edge2=-Inf)

windows()
ggplot(AUC.long, aes(species, value, group=type, fill=type))+geom_bar(stat="identity", position="dodge")+theme_classic()+geom_rect(data=thresholds, aes(xmin=edge1, xmax=edge2, ymin=0, ymax=inaccurate), fill="#d8b365", colour="brown",alpha=.2, inherit.aes = FALSE)+geom_rect(data=thresholds, aes(xmin=edge1, xmax=edge2, ymin=inaccurate, ymax=accurate), fill="#5ab4ac", colour="#01665e",alpha=.2, inherit.aes = FALSE)+geom_rect(data=thresholds, aes(xmin=edge1, xmax=edge2, ymin=accurate, ymax=1), fill="#5ab4ac", alpha=.5, colour="#01665e",inherit.aes = FALSE)+geom_bar(colour="black",stat="identity", position="dodge")+theme_classic()+coord_flip()+ylab("AUC")+xlab("")+theme(legend.position = "bottom")+scale_fill_manual(name="", values=c("grey", "white"))
ggsave("AUC_conditional_marginal_streams.png", height=12, width=8, units="in")
