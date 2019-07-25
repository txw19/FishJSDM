# rm(list=ls())
library(ggplot2)
library(reshape2)
library(ggrepel)

setwd("~/MN fish trends/Joint Distribution Model/JDM_Fish_Intro_paper/Model_devel_gjam_lakes/")
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



 sppnames <- c("black bullhead","black crappie","bluegill", "bowfin", "brown bullhead", "common carp","golden shiner","green sunfish","hybrid sunfish","largemouth bass","northern pike","pumpkinseed", "rock bass","smallmouth bass","walleye","white sucker","yellow bullhead","yellow perch","coldwater")


colnames(ERhoMean)=sppnames
rownames(ERhoMean)=sppnames 

corr.data=melt(ERhoMean, value.name = "Env.cor", varnames = c("Species1", "Species2"))
 
colnames(RhoMean)=sppnames
rownames(RhoMean)=sppnames  
 
corr.data2=melt(RhoMean, value.name = "Res.cor", varnames = c("Species1", "Species2"))

corr.data=merge(corr.data, corr.data2, by=c("Species1", "Species2"))
corr.data=corr.data[complete.cases(corr.data),]
#if going to facet, need to repeat so all species dots appear in all species panels

corr.data.doubled=melt(ERhoMean, value.name = "Env.cor", varnames = c("Species2", "Species1"))
corr.data.doubled2=melt(RhoMean, value.name = "Res.cor", varnames = c("Species2", "Species1"))
corr.data.doubled=merge(corr.data.doubled, corr.data.doubled2, by=c("Species1", "Species2"))
corr.data.doubled=corr.data.doubled[complete.cases(corr.data.doubled),]
corr.data.doubled=rbind(corr.data.doubled, corr.data)

windows()
ggplot(corr.data.doubled, aes(Env.cor, Res.cor, colour=Species1))+geom_hline(yintercept = 0, alpha=0.5, colour="grey")+geom_vline(xintercept=0, alpha=0.5, colour="grey")+geom_point()+geom_text_repel(aes(label=Species1), size=3) +facet_wrap(~Species2)+theme_bw()+xlab("Shared environmental correlation")+ylab("Residual correlation")+theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank())
ggsave("environmental_and_residual_correlations_by_species_PA.png", height=12, width=10, units="in")

 