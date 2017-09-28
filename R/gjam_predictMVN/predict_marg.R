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


# !!!!!!!!!!! TO DO !!!!!!!!!!!
# Standardize predictors (gjam does this internally, we need to do so here?)



# Read in RDS files in order, 1-10
out1 <- readRDS(file="gjamOut1.rds")
out2 <- readRDS(file="gjamOut2.rds")
out3 <- readRDS(file="gjamOut3.rds")
out4 <- readRDS(file="gjamOut4.rds")
out5 <- readRDS(file="gjamOut5.rds")
out6 <- readRDS(file="gjamOut6.rds")
out7 <- readRDS(file="gjamOut7.rds")
out8 <- readRDS(file="gjamOut8.rds")
out9 <- readRDS(file="gjamOut9.rds")
out10 <- readRDS(file="gjamOut10.rds")
# Put gjam MCMC output into list
mcmc_list <- list(out1,out2,out3,out4,out5,out6,out7,out8,out9,out10)


########## GET BETA MCMC SAMPLES INTO CORRECT FORMAT
## !!!!!!!!PLACE THIS INSIDE FOR LOOP ONCE WORKING
# Beta coefficients dim = [nsim, 96], need array of [nsim, nspecies, ncovariates]
Beta <- out1$chains$bgibbs #  will be list[[i]]$chains$bgibbs below
dim(Beta)

# Remove intercepts from Beta

int.names <- c("black.bullhead_intercept","black.crappie_intercept","bowfin.(dogfish)_intercept","brown.bullhead_intercept",
               "common.carp_intercept","golden.shiner_intercept","largemouth.bass_intercept","northern.pike_intercept",
               "rock.bass_intercept","smallmouth.bass_intercept","tullibee.(cisco)_intercept","walleye_intercept",
               "white.sucker_intercept","yellow.bullhead_intercept","yellow.perch_intercept","sunfish_intercept")

Beta2 <- Beta[, !colnames(Beta) %in% int.names]

dim(Beta2)
head(Beta2)

# We need matrices, one for each covariate, with rows=num MCMC samples, columns = number of species
# and a list for each covariates. 
const <- 5
lst1 <- list()
for(i in 1:5){
  lst1[[i]] <- Beta2sub[,c(i,5+i,10+i,15+i,20+i,25+i,30+i,35+i,40+i,45+i,50+i,55+i,60+i,65+i,70+i,75+i)]
}
str(lst1)
head(lst1[[1]])

BetaA <- array(as.numeric(unlist(lst1)), dim=c(length(ID),numY,K))

dim(BetaA)

######### GET SIGMA IN CORRECT DIMENSIONS - PLACE IN FOR LOOP ONCE WORKIN
# PUT IN LIST AS ABOCE AND THEN INTO AN ARRAY OF APPROPRIATE DIMENSIONS
# SigmaOut.rds contains the lower triangle of the variance covariance matrix, i.e., 136 columns with 16 species
Sigma <- out1$chains$sgibbs
str(Sigma)
dim(Sigma)

numY <- 4

postSigMeans <- apply(Sigma,2,mean)
length(postSigMeans)

# S <- diag(numY)
# S[lower.tri(S, diag=TRUE)] <- postSigMeans

# Convert S (half matrix) to full matrix
X <- diag(numY)
X[lower.tri(X, diag=TRUE)] <- postSigMeans
X <- X + t(X) - diag(diag(X))

# Convert Sigma to Rho
K <- numY # Number of species
Rho <- array(NA, dim=c(dim(Sigma)[1], dim(X)[1], dim(X)[1]))
dim(Rho)
for(samp in 1:dim(Sigma)[1]){
  # Grab mcmc sample from Sigma
  sigsamp <- Sigma[samp,]
  # Convert S (half matrix) to full matrix
  Stemp <- diag(16) 
  Stemp[lower.tri(Stemp, diag=TRUE)] <- sigsamp 
  Stemp <- Stemp + t(Stemp) - diag(diag(Stemp)) 
  # Calculate correlation matrix from Sigma
  for (k in 1:K){
    for (k.prime in 1:K){
      Rho[samp, k, k.prime] <- Stemp[k, k.prime]/
        sqrt(Stemp[k,k]*Stemp[k.prime,k.prime])
    }
  }  
}














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
length(ID)

## Loop over all 10 MCMC outs and predict out-of-sample dataset to calculate RMSE
# Container to hold predicted values
predict_list <- list()
# List for holding marginal posterior predictions
marg_pred_posterior <- list()
###### Start loop
for(m in 1:length(mcmc_list)){

####### Grab posterior samples for all parameters
# Intercept 
a <- mcmc_list[[m]]$alpha[ID,]
# Slopes
for(i in 1:dim(mcmc_list[[1]]$beta)[3]){
  assign(paste("b", i, sep = ""), mcmc_list[[m]]$beta[ID, ,i])
}
#######!!!!!!!!! Requires manual placing of b1, b2, ... bn at the moment
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
datOOS <- datFinal[index == m,] 
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

# ypred dim := # iterations, # responses, # lakes to predict
ypred <- array(NA, c(dim(a)[1] , numY, dim(newX)[1]))
for(i in 1:dim(a)[1]){ # iterations
  sigma_mat <- matrix(unlist(mcmc_list[[m]]$Sigma.B[i,,]), ncol = numY, byrow = TRUE) # Grab mcmc samples of sigma and turn them into matrix
  for(k in 1:dim(newX)[1]){ # lakes to predict
    ypred[i,,k] <- rmvnorm(1, newX[k,]%*%betam[i,,], sigma_mat)
  }
}

# Hold all MCMC posterior predicted values for all lakes
marg_pred_posterior[[m]] <- ypred

# Marginal posterior predicted mean for each lake and response variable
marg.predictions <- apply(ypred,c(3,2),mean)

pred_obs <- cbind(marg.predictions, newY)

# Save predicted values and observed
predict_list[[m]] <- pred_obs

} #####!!!!!!!!!! End mcmc_list for loop

# # Marginal posterior predicted mean for each lake and response variable
# 
# marg.predictions <- apply(ypred,c(3,2),mean)
# # write.csv(marg.predictions,"marg.predictions.csv")

# Calculate RMSE
# Y1 - Y5 (log_tp, log_tn, log_chla, log_secchi, log_no3)

head(predict_list[[1]])
# Converve predict_list into matrix
predict_matrix <- do.call(rbind, predict_list)
  # matrix(unlist(predict_list), ncol = 10, byrow = F)
dim(predict_matrix)
# head(predict_list[[1]])
# head(predict_matrix)
# tail(predict_list[[10]])
# tail(predict_matrix)

# Seperate predictions and observed for rmse calcs
marg.predictions <- predict_matrix[,1:5]
obs.Y <- predict_matrix[,6:10]

rmse <- function(observed, predicted){
  sqrt( sum( (observed - predicted)^2 , na.rm = TRUE ) / sum( !is.na(observed) ) )
}

marg_rmse_calcs <- matrix(NA, ncol=5, nrow=1)
for(i in 1:numY){
  marg_rmse_calcs[i] <- rmse(obs.Y[,i],marg.predictions[,i])
}  
marg_rmse_calcs

colnames(marg_rmse_calcs) <- c('log(TP)', 'log(TN)','log(CHL)','log(Secchi)','log(NO3')

write.csv(marg_rmse_calcs, 'marginal_rmse.csv', row.names = F)


# Calculate R2
r2_calc <- function(observed, predicted){
  1-(sum( (observed - predicted)^2 , na.rm = TRUE ) / sum( (observed - mean(observed, na.rm=T))^2 , na.rm = TRUE ))
}

r2_calcs <- matrix(NA, ncol=5, nrow=1)
for(i in 1:numY){
  r2_calcs[i] <- r2_calc(obs.Y[,i],marg.predictions[,i])
}  
r2_calcs

colnames(r2_calcs) <- c('log(TP)', 'log(TN)','log(CHL)','log(Secchi)','log(NO3')
write.csv(r2_calcs, 'marginal_r2.csv', row.names = F)



# Save predicted posterior distibutions
saveRDS(marg_pred_posterior, file="marg_post_pred.rds")

##############  


# # Erin code
# predict.all=function(run,b,e,iters,newX){
#   Y.pred=array(dim=c(nrow(newX),ncol(run$data$Y),iters))
#   s=sample(b:e,iters)
#   for(j in 1:iters){
#     for(i in 1:nrow(newX)){
#       Y.pred[i,,j]=rmvnorm(1,newX[i,]%*%t(run$beta[,,s[j]]), run$Sigma[,,s[j]])	
#     }
#   }
#   out=list()
#   out$Y.pred=Y.pred
#   out$newX=newX
#   return(out)
# }
