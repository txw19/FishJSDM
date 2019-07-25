# rm(list=ls())
library(dplyr)
library(ggplot2)


# Read in saved data and compute stats
betaout <- readRDS(file="BetaOut.rds")

finalDat=read.csv("finalDat.csv")

sppnames <- colnames(finalDat)[8:24]
sppnames[17]="coldwaterspp"

remove.names=  paste( sppnames,"intercept", sep = "_")
# Remove intercepts from Beta

betaout <- betaout[, !colnames(betaout) %in% remove.names]
dim(betaout)

#make dataframe for ggplot

betaMeans <- apply(betaout,2,mean)

effects=data.frame("effects"=names(betaMeans))
toplot = transform(effects, variable = colsplit(effects, pattern = "_", names = c('1', '2')))

toplot$mean.beta=as.numeric(betaMeans)


toplot$betaUCI <- apply(betaout,2,quantile,  0.975)
toplot$betaLCI <- apply(betaout,2,quantile,  0.025)

toplot$significant="N"
toplot$significant[toplot$betaUCI*toplot$betaLCI>0]="Y"


#####################################################
########### PLOT ####################################
#####################################################

windows()
ggplot(toplot, aes(variable.2, mean.beta, colour=significant, shape=significant))+geom_pointrange(aes(ymax=betaUCI, ymin=betaLCI))+facet_wrap(~variable.1)+coord_flip()+theme_bw()+geom_hline(yintercept=0)+scale_colour_manual(values=c("grey", "royalblue"))+scale_shape_manual(values=c(1,19))+theme(legend.position="none",strip.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("Estimated effect")+xlab("")

ggsave("lakes_effect_estimates_by_species.png", height=8, width=12, units="in")


#facet by effect instead
windows()
ggplot(toplot, aes(variable.1, mean.beta, colour=significant, shape=significant))+geom_pointrange(aes(ymax=betaUCI, ymin=betaLCI))+facet_wrap(~variable.2)+coord_flip()+theme_bw()+geom_hline(yintercept=0)+scale_colour_manual(values=c("grey", "royalblue"))+scale_shape_manual(values=c(1,19))+theme(legend.position = "none", strip.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("Estimated effect")+xlab("")

ggsave("lakes_effect_estimates_by_effect.png", height=8, width=8, units="in")

