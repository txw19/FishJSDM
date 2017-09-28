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


###### Create the out-of-sample test data set #######
# Index for splitting dat into 10 oos datasets (~10% of the total), will allow each lake will get predicted oos
ind <- rep(1:10,dim(dat)[1]/10)
ind <- c(ind, c(1,2,3))
length(ind)
# Randomize index
ind <- sample(ind)

dat$index <- ind
table(dat$index)
summary(dat)
# Write out final dataset with OOS index appended
write.csv(dat, "dat.csv", row.names = F)
dim(dat)

########### Begin loop through datasets for 10-fold cross-validation
# Container for predictions
preds <- list()
# Loop through oos datasets
for(i in 1:10){
  
  # # Select out dataset for model fitting (dat_gjam)
  # dat_gjam <- dat[index != 1,]
  # dim(dat_gjam)


xdat <- dat[,20:24]
# sqrt-transform CPE
ydat <- sqrt(dat[,c(2,4,8,17)]) # Grab a few spp for testing
head(ydat)

# Some gjam options for modelList:
# holdoutN = 0, number of observations to hold out for out-of-sample prediction.
# holdoutIndex =  numeric(0), numeric vector of observations (row numbers) to holdout for out-of-sample prediction

# Rows for holdoutIndex for 10-fold x-validation (these are used in gjam for marginal prediction for oos data??)
oos_rows <- which(dat$index==i)

# -typeNames = "CA" means continuous abundance (meaning greater than 0).
start.time = Sys.time()  # Start timer 

jdm1 = gjam(~ log_area + log_depth + Secchi.lake.mean + mean.gdd + alkalinity, xdata=xdat, ydata=ydat,
          modelList=list(PREDICTX = F, ng=50000,burnin=20000,typeNames=rep("CA",4), holdoutIndex=oos_rows))

end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

# summary(jdm1)
# gjamPredict(jdm1, y2plot = colnames(ydat)) 

#####################################################
##### MARGINAL predictions for oos data
#####################################################
# Predicted values of Y
yMu  <- jdm1$prediction$ypredMu    
#holdout observations (rows)
hold <- jdm1$modelList$holdoutIndex    

# Observed Y for oos samples
ObsY <- jdm1$inputs$y[hold,]
# Predicted Y for oos samples
MargPredY <- yMu[hold,]

# Plot observed vs. predicted
# plot(ObsY, PredY, xlab='True', ylab='')

######################################
# out-of-sample (MARGINAL) prediction can also be done usig gjamPredict
# xdata <- xdat[hold,]
# newdata   <- list(xdata = xdata, nsim = 500 )
# p1 <- gjamPredict(jdm1, newdata = newdata)
# # str(p1)
# 
# plot(ObsY, p1$sdList$yMu)

############################################
### CONDITIONAL predictions for oos data
############################################

ydataCond <- jdm1$inputs$y[hold,]
# dim(ydataCond)
# head(ydataCond)
xdata <- xdat[hold,]

newdata <- list(xdata = xdata, condY = ydataCond, nsim=2000)

condPreds <- gjamPredict(output = jdm1, newdata = newdata)
# str(condPreds)
CondPredY <- condPreds$sdList$yMu

# Plot observed vs. predicted
# plot(ObsY, CondPredY, xlab='True', ylab='')
##############################################
# Compile marginal and conditional predictions
##############################################
# Put observed and predicted in list
preds[[i]] <- cbind(MargPredY,CondPredY,ObsY)

# preds <- cbind(MargPredY,CondPredY,ObsY)

## Export all 10 model fit summaries (beta's, sigma, and gjam fit objects)
# these are all the betas
BetaOut <- jdm1$chains$bgibbs
# sigma
SigmaOut <-jdm1$chains$sgibbs
saveRDS(BetaOut, file=paste("BetaOut",i,".rds",sep=''))
saveRDS(SigmaOut, file=paste("SigmaOut",i,".rds",sep=''))
saveRDS(jdm1, file=paste("gjamOut",i,".rds",sep=''))

} # end for loop



## NOT RUN - now deriving predictived values using custom script
####################################
## Calculate RMSPE and predictive R2
####################################
# Number of species
# numY <- 16
# 
# Converve preds - a list - into a matrix
# Export gjam predictions for comparison with ours
pred_matrix <- do.call(rbind, preds)
dim(pred_matrix)

write.csv(gjam_pred_matrix, "PosteriorPreds.csv", row.names = F)
# 
# # Seperate predictions and observed for rmse calcs
# marg.predictions <- pred_matrix[,1:16]
# cond.predictions <- pred_matrix[,17:32]
# obs.Y <- pred_matrix[,33:48]
# 
# # rmse function
# rmse <- function(observed, predicted){
#   sqrt( sum( (observed - predicted)^2 , na.rm = TRUE ) / sum( !is.na(observed) ) )
# }
# 
# 
# # conditional rmspe
# cond_rmse_calcs <- matrix(NA, ncol=5, nrow=1)
# for(i in 1:numY){
#   cond_rmse_calcs[i] <- rmse(obs.Y[,i],cond.predictions[,i])
# }  
# cond_rmse_calcs
# 
# cond_rmse_calcs <- matrix(cond_rmse_calcs, nrow=1, ncol=numY)
# 
# colnames(cond_rmse_calcs) <- c("black.bullhead_intercept","black.crappie_intercept","bowfin.(dogfish)_intercept","brown.bullhead_intercept",
#                                             "common.carp_intercept","golden.shiner_intercept","largemouth.bass_intercept","northern.pike_intercept",
#                                             "rock.bass_intercept","smallmouth.bass_intercept","tullibee.(cisco)_intercept","walleye_intercept",
#                                             "white.sucker_intercept","yellow.bullhead_intercept","yellow.perch_intercept","sunfish_intercept")
# 
# write.csv(cond_rmse_calcs, 'conditional_rmse.csv', row.names = F)
# 
# 
# # Marginal rmspe
# marg_rmse_calcs <- matrix(NA, ncol=5, nrow=1)
# for(i in 1:numY){
#   marg_rmse_calcs[i] <- rmse(obs.Y[,i],marg.predictions[,i])
# }  
# marg_rmse_calcs
# 
# marg_rmse_calcs <- matrix(marg_rmse_calcs, nrow=1, ncol=numY)
# 
# colnames(marg_rmse_calcs) <- c("black.bullhead_intercept","black.crappie_intercept","bowfin.(dogfish)_intercept","brown.bullhead_intercept",
#                               "common.carp_intercept","golden.shiner_intercept","largemouth.bass_intercept","northern.pike_intercept",
#                               "rock.bass_intercept","smallmouth.bass_intercept","tullibee.(cisco)_intercept","walleye_intercept",
#                               "white.sucker_intercept","yellow.bullhead_intercept","yellow.perch_intercept","sunfish_intercept")
# 
# write.csv(marg_rmse_calcs, 'marginal_rmse.csv', row.names = F)
# 
# ###############################################
# # Calculate conditional predictive R2
# r2_calc <- function(observed, predicted){
#   1-(sum( (observed - predicted)^2 , na.rm = TRUE ) / sum( (observed - mean(observed, na.rm=T))^2 , na.rm = TRUE ))
# }
# 
# r2_calcs <- matrix(NA, ncol=5, nrow=1)
# for(i in 1:numY){
#   r2_calcs[i] <- r2_calc(obs.Y[,i],cond.predictions[,i])
# }  
# r2_calcs
# 
# r2_calcs <- matrix(r2_calcs, nrow=1, ncol=numY)
# 
# colnames(r2_calcs) <- c("black.bullhead_intercept","black.crappie_intercept","bowfin.(dogfish)_intercept","brown.bullhead_intercept",
#                         "common.carp_intercept","golden.shiner_intercept","largemouth.bass_intercept","northern.pike_intercept",
#                         "rock.bass_intercept","smallmouth.bass_intercept","tullibee.(cisco)_intercept","walleye_intercept",
#                         "white.sucker_intercept","yellow.bullhead_intercept","yellow.perch_intercept","sunfish_intercept")
# 
# write.csv(r2_calcs, 'cond_r2.csv', row.names = F)
# 
# 
# ##################################
# # Calculate marginal predictive R2
# r2_calcs_marg <- matrix(NA, ncol=5, nrow=1)
# for(i in 1:numY){
#   r2_calcs_marg[i] <- r2_calc(obs.Y[,i],marg.predictions[,i])
# }  
# r2_calcs_marg
# 
# r2_calcs_marg <- matrix(r2_calcs_marg, nrow=1, ncol=numY)
# 
# colnames(r2_calcs_marg) <- c("black.bullhead_intercept","black.crappie_intercept","bowfin.(dogfish)_intercept","brown.bullhead_intercept",
#                         "common.carp_intercept","golden.shiner_intercept","largemouth.bass_intercept","northern.pike_intercept",
#                         "rock.bass_intercept","smallmouth.bass_intercept","tullibee.(cisco)_intercept","walleye_intercept",
#                         "white.sucker_intercept","yellow.bullhead_intercept","yellow.perch_intercept","sunfish_intercept")
# 
# write.csv(r2_calcs_marg, 'marginal_r2.csv', row.names = F)
# 
# 
