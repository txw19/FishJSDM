# rm(list=ls())
require(data.table)
require(mvtnorm)

### USER INPUT
# 1. Change nsim to dictate the number of MCMC samples used for prediction
# 2. Fill in the betam arrary with the appropriate number of b's corresponding to predictors
# 3. Create the appropriate newX matrix to correspond to the predictors used when fitting the model

# Predict Y1 - Yn

# Read in final dat
datFinal <- fread("dat.csv")

dim(datFinal)
# Read in RDS files in order, 1-10
out1 <- readRDS(file="JDMmcmc_out1.rds")
out2 <- readRDS(file="JDMmcmc_out2.rds")
out3 <- readRDS(file="JDMmcmc_out3.rds")
out4 <- readRDS(file="JDMmcmc_out4.rds")
out5 <- readRDS(file="JDMmcmc_out5.rds")
out6 <- readRDS(file="JDMmcmc_out6.rds")
out7 <- readRDS(file="JDMmcmc_out7.rds")
out8 <- readRDS(file="JDMmcmc_out8.rds")
out9 <- readRDS(file="JDMmcmc_out9.rds")
out10 <- readRDS(file="JDMmcmc_out10.rds")
# Put MCMC output into list
mcmc_list <- list(out1,out2,out3,out4,out5,out6,out7,out8,out9,out10)

# names(mcmc_list[[1]])

# Number of response variables
numY <- dim(mcmc_list[[1]]$alpha)[2]
####!!! Number of desired samples for prediction
nsim <- 1000
# Chain length from analysis
chainLength <- dim(mcmc_list[[1]]$alpha)[1]
# # Select thinned steps in chain for posterior predictions to ensure we take values from 
# # length of posterior
ID = seq( 1 , chainLength , floor(chainLength/nsim) )

## Loop over all 10 MCMC outs and predict out-of-sample dataset to calculate RMSE
# Container to hold predicted values
predict_cond_list <- list()
# List for holding conditional posterior predictions
cond_pred_posterior <- list()

###### Start loop
for(m in 1:length(mcmc_list)){

####### Grab posterior samples for all parameters
# Intercept 
a <- mcmc_list[[m]]$alpha[ID,]
# Slopes
for(i in 1:dim(mcmc_list[[1]]$beta)[3]){
  assign(paste("b", i, sep = ""), mcmc_list[[m]]$beta[ID, ,i])
}
####### Requires manual placing of b1, b2, ... bn at the moment
# Take mcmc samples of intercept and slopes and put them in array for matrix multiplication
# Number of parameters (slopes and intercept (+1)) estimated
numParams <- (dim(mcmc_list[[1]]$beta)[3] + 1)
## Create new array of coefficients for matrix multiplication
betam <- array(NA, c(dim(a)[1],numParams,numY )) # numY = the number of response variables
for(i in 1:numY){
  betam[,,i] <- c(a[,i], b1[,i], b2[,i], b3[,i],b4[,i],b5[,i],b6[,i],b7[,i],b8[,i])
}

# Covaraince matrix
tsigma <- mcmc_list[[m]]$Sigma.B[ID,,]
#########

# Grab oos dataset for new X's
datOOS <- datFinal[index ==m,] 
####### !!!!!!!! Need to specify predictors
# newX <- cbind(datOOS$log_depth, datOOS$log_area)
newX <- as.matrix(datOOS[,71:78])
# Add 1's to design matrix for matrix multiplication
newX <- as.matrix(cbind(c(rep(1,dim(newX)[1])),newX))
# Grab new Y's
newY <- datOOS[, c('log_tp','log_tn','log_chla','log_secchi','log_no3'), with=FALSE]
# Convert newY to a matrix for below calculations
newY <- as.matrix(newY)
colnames(newY) <- NULL

# Y.cond=array(dim=c(nrow(newX),ncol(run$data$Y),iters))
# newx, numY, iters
# Erin y.cond[i,,j]
# y.cond[newX, numY, iters]

# Ty
# y.cond[iters, numY, newX]

Y.cond=array(NA, c(dim(a)[1] , numY, dim(newX)[1]))
for(k in 1:nrow(newX)){
  f=which(is.na(newY[k,])!=T) # Which nutrients have data? Which are NOT missing?
  if(length(f)==0){ # if all Y's are missing predict Y.cond as follows, i.e., there is nothing to condition on
    for(i in 1:dim(a)[1]){
      sigma_mat <- matrix(unlist(mcmc_list[[m]]$Sigma.B[i,,]), ncol = numY, byrow = TRUE) # Grab mcmc samples of sigma and turn them into matrix
      Y.cond[i,,k]=rmvnorm(1, newX[k,]%*%betam[i,,], sigma_mat)}}
  if(length(f)!=0){ # if there are some observed nutrients, then...
    for(i in 1:dim(a)[1]){
      C=matrix(unlist(mcmc_list[[m]]$Sigma.B[i,,]), ncol = numY, byrow = TRUE) # Grab Sigma matrix sample
      Mu=newX[k,]%*%betam[i,,] # Calculte parameter vector Mu of MVN
      if(length(f)!=numY){ # Predict unobserved Y's, given some others observed - if the number of nutrients observed is not 4 (the total number of nutrients in this example), then 
        Y.cond[i,-f,k]=rmvnorm(1,mean=Mu[-f]+C[-f,f]%*%solve(C[f,f])%*%(newY[k,f]-Mu[f]), 
                               sigma=C[-f,-f]-C[-f,f]%*%solve(C[f,f])%*%C[f,-f])}
      if(length(f)==1){Y.cond[i,f[1],k]=rnorm(1,Mu[f[1]],sqrt(C[f[1],f[1]]))} # If there is a single nutrient then there is nothing to condition on, predict based on marginal variance 
      if(length(f)>1){
        for(fs in 1:length(f)){ # Predict observed Y's, conditional on others observed
          Y.cond[i,f[fs],k]=rmvnorm(1,
                                    mean=Mu[f[fs]]+C[f[fs],f[-fs]]%*%solve(C[f[-fs],f[-fs]])%*%(newY[k,f[-fs]]-Mu[f[-fs]]),
                                    sigma=C[f[fs],f[fs]]-C[f[fs],f[-fs]]%*%solve(C[f[-fs],f[-fs]])%*%C[f[-fs],f[fs]])
        }
      }
    }}}

# Hold all MCMC posterior predicted values for all lakes
cond_pred_posterior[[m]] <- Y.cond

# Posterior predicted mean for each lake and response variable
cond.predictions <- apply(Y.cond,c(3,2),mean)
# Put predictions and observations together
pred_cond_obs <- cbind(cond.predictions, newY)
# Store predictions and obs in list
predict_cond_list[[m]] <- pred_cond_obs

} ##### End mcmc_list for loop

# Calculate RMSE
# Y1 - Y5 (log_tp, log_tn, log_chla, log_secchi, log_no3)
head(predict_cond_list[[1]])
# Converve predict_list into matrix
predict_cond_matrix <- do.call(rbind, predict_cond_list)
# matrix(unlist(predict_list), ncol = 10, byrow = F)
dim(predict_cond_matrix)

# Seperate predictions and observed for rmse calcs
cond.predictions <- predict_cond_matrix[,1:5]
obs.Y <- predict_cond_matrix[,6:10]

rmse <- function(observed, predicted){
  sqrt( sum( (observed - predicted)^2 , na.rm = TRUE ) / sum( !is.na(observed) ) )
}

cond_rmse_calcs <- matrix(NA, ncol=5, nrow=1)
for(i in 1:numY){
  cond_rmse_calcs[i] <- rmse(obs.Y[,i],cond.predictions[,i])
}  
cond_rmse_calcs

colnames(cond_rmse_calcs) <- c('log(TP)', 'log(TN)','log(CHL)','log(Secchi)','log(NO3')

write.csv(cond_rmse_calcs, 'conditional_rmse.csv', row.names = F)


# Calculate R2
r2_calc <- function(observed, predicted){
  1-(sum( (observed - predicted)^2 , na.rm = TRUE ) / sum( (observed - mean(observed, na.rm=T))^2 , na.rm = TRUE ))
}

r2_calcs <- matrix(NA, ncol=5, nrow=1)
for(i in 1:numY){
  r2_calcs[i] <- r2_calc(obs.Y[,i],cond.predictions[,i])
}  
r2_calcs

colnames(r2_calcs) <- c('log(TP)', 'log(TN)','log(CHL)','log(Secchi)','log(NO3')
write.csv(r2_calcs, 'cond_r2.csv', row.names = F)




# Save predicted posterior distibutions
saveRDS(cond_pred_posterior, file="cond_post_pred.rds")


#################################
################# TEST CRAP
#################################
# head(newY)
# f=which(is.na(newY[4,])!=T)  # obs 1, don't observe TN (element 2)
# f
# sigma_mat <- matrix(unlist(mcmc_list[[1]]$Sigma.B[1,,]), ncol = numY, byrow = TRUE)
# # C[f[fs],f[-fs]]%*%solve(C[f[-fs],f[-fs]])
# # C[f[fs],f[-fs]]=
# sigma_mat[f[1],f[-1]]
# # C[f[-fs],f[-fs]] =
# sigma_mat[f[-1],f[-1]]

# for(fs in 1:length(f)){
#   Y.cond[i,f[fs],j]=rmvnorm(1,
#                             mean=Mu[f[fs]]+C[f[fs],f[-fs]]%*%solve(C[f[-fs],f[-fs]])%*%(run$data$Ynew[i,f[-fs]]-Mu[f[-fs]]),
#                             sigma=C[f[fs],f[fs]]-C[f[fs],f[-fs]]%*%solve(C[f[-fs],f[-fs]])%*%C[f[-fs],f[fs]])
# }
# }
# log_tp   log_tn log_chla  log_secchi  log_no3
# 2.014903    NA   1.360977  1.80828877 2.995732

# f
# 
# mu <- c(0.5, 0.6, 0.7, 0.8, 0.9)
# 
# 
# # if(length(f)!=5)
# mu[-f] # (unobserved) mu
# mu[f] # Observed mu's
# sigma_mat[-f,f]
# sigma_mat[f,f]
# 
# 
# # if(length(f)>1)
# fs=2
# mu[f[fs]] # observed
# mu[f[-fs]]



##################################
# Erin's code
# predict.cond=function(run,b,e,iters,newX){
#   Y.cond=array(dim=c(nrow(newX),ncol(run$data$Y),iters))
#   s=sample(b:e,iters) # Take every nth mcmc sample
#   for(i in 1:nrow(newX)){
#     f=which(is.na(run$data$Ynew[i,])!=T) # Which nutrients have data? Which are NOT missing
#     if(length(f)==0){ # if all Y's are missing predict Y.cond as follows, i.e., there is nothing to condition on
#       for(j in 1:iters){
#         Y.cond[i,,j]=rmvnorm(1,newX[i,]%*%t(run$beta[,,s[j]]),run$Sigma[,,s[j]])}}
#     if(length(f)!=0){ # if there are some observed nutrients, then...
#       for(j in 1:iters){
#         C=run$Sigma[,,s[j]] # Grab Sigma matrix sample
#         Mu=newX[i,]%*%t(run$beta[,,s[j]]) # Calculte parameter vector Mu of MVN
#         if(length(f)!=4){ # Predict unobserved Y's, given some others observed - if the number of nutrients observed is not 4 (the total number of nutrients in this example), then 
#           Y.cond[i,-f,j]=rmvnorm(1,mean=Mu[-f]+C[-f,f]%*%solve(C[f,f])%*%(run$data$Ynew[i,f]-Mu[f]), 
#                                  sigma=C[-f,-f]-C[-f,f]%*%solve(C[f,f])%*%C[f,-f])}
#         if(length(f)==1){Y.cond[i,f[1],j]=rnorm(1,Mu[f[1]],sqrt(C[f[1],f[1]]))} # If there is a single nutrient then there is nothing to condition on, predict based on marginal variance 
#         if(length(f)>1){
#           for(fs in 1:length(f)){ # Predict observed Y's, conditional on others observed
#             Y.cond[i,f[fs],j]=rmvnorm(1,
#                                       mean=Mu[f[fs]]+C[f[fs],f[-fs]]%*%solve(C[f[-fs],f[-fs]])%*%(run$data$Ynew[i,f[-fs]]-Mu[f[-fs]]),
#                                       sigma=C[f[fs],f[fs]]-C[f[fs],f[-fs]]%*%solve(C[f[-fs],f[-fs]])%*%C[f[-fs],f[fs]])
#           }
#         }
#       }}}
# }
# 





