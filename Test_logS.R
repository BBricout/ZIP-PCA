rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

setwd("~/Documents/ZIP-PCA") # A changer
library(Rcpp)
library(PLNmodels)

source("codesBB/UtilsBB.R")
source("Data_create.R")
sourceCpp("src/optim_rank_zip_log_new.cpp")
sourceCpp("src/optim_rank_ZIP_new.cpp")

#########################################################
######## Paramètres initiaux ###########################
########################################################


init <- InitZiPLN(data,q) 
params <- list(B = as.matrix(init$mStep$beta), 
               D = as.matrix(init$mStep$gamma), 
               C = as.matrix(init$mStep$C), 
               M = as.matrix(init$eStep$M), 
               S = as.matrix(init$eStep$S))

paramsLogS <- list(B = as.matrix(init$mStep$beta), 
               D = as.matrix(init$mStep$gamma), 
               C = as.matrix(init$mStep$C), 
               M = as.matrix(init$eStep$M), 
               logS = log(as.matrix(init$eStep$S)))
config <- PLNPCA_param()$config_optim ; tolXi <- 1e-04

parmsLogSInit <- c(init$mStep$gamma, init$mStep$beta, as.vector(init$mStep$C), 
                   as.vector(t(init$eStep$M)), log(as.vector(t(init$eStep$S))))


################################################################################
## Tests sur les calculs d'Elbo et de gradient


obj_ref <- Elbo_grad_Rcpp(data, params, tolXi)
obj_log <- Elbo_grad_logS_Rcpp(data, paramsLogS, tolXi)

BobjLogS <- function(parmsLogS, data){
  params <- Parms2Params(parms=parmsLogS, data=data)
  params$S <- exp(params$S)
  Elbo_grad_Rcpp(data=data, params=params, tolXi)$obj
}
BgradLogS <- function(parmsLogS, data){
  params <- Parms2Params(parms=parmsLogS, data=data)
  params$S <- exp(params$S)
  elbo <- Elbo_grad_Rcpp(data=data, params=params, tolXi)
  c(as.vector(elbo$gradD), as.vector(elbo$gradB), as.vector(elbo$gradC),
    as.vector(t(elbo$gradM)), as.vector(t(elbo$gradS))*as.vector(t(params$S)))
}

obj_Ste <- BobjLogS(parmsLogSInit, data)

grad_Ste <- BgradLogS(parmsLogSInit, data)
grad_ref <- Parms2Params(grad_Ste, data)

plot(obj_ref$gradB, obj_log$gradB) ; abline(0,1)
plot(obj_ref$gradC, obj_log$gradC) ; abline(0,1)
plot(obj_ref$gradD, obj_log$gradD) ; abline(0,1)
plot(obj_ref$gradM, obj_log$gradM) ; abline(0,1)
plot(grad_ref$S, obj_log$gradS) ; abline(0,1)

################################################################################
############### Tests de l'algo ##############################################
###########################################################################

out <- nlopt_optimize_ZIP_logS(data, paramsLogS, config, tolXi)


#### Estimations

B.hat <- out$B
D.hat <- out$D
C.hat <- out$C
M.hat <- out$M
S.hat <- out$S

XB <- VectorToMatrix(data$X %*% true$mStep$beta, n, p)
XB.hat <- VectorToMatrix(data$X %*% B.hat, n, p)
XB.init <- VectorToMatrix(data$X %*% params$B, n, p)


XD <- VectorToMatrix(data$X %*% true$mStep$gamma, n, p)
XD.hat <- VectorToMatrix(data$X %*% D.hat, n, p)
XD.init <- VectorToMatrix(data$X %*% params$D, n, p)


pred.init <- pred <- exp(XB.init) #  + params$M %*% t(params$C) + 0.5 * (params$S*params$S) %*% t(params$C * params$C))
pred <- exp(XB.hat + M.hat %*% t(C.hat) + 0.5 * (exp(S.hat)*exp(S.hat)) %*% t(C.hat * C.hat))

Y.logit <- ifelse(data$Y == 0, 0, 1)

############################ Plots

#---------------------Elbo

plot(out$objective_values[out$objective_values < 1*out$objective[1]],type='l', main = "ELBO path",ylab='ELBO',xlab='iter')

#----------------- Estimations

par(mfrow=c(1, 2))
plot(XD, XD.init, col='red', main = "Estimation de la logistique") ; abline(0,1)
points(XD, XD.hat) ; abline(0,1)



plot(XB, XB.init,col='red', main = "Estimation des régresseurs") ; abline(0,1)
points(XB, XB.hat) ; abline(0,1)


par(mfrow=c(1,3))
boxplot(VectorToMatrix(XD, n, p) ~ Y.logit, main = "True")
boxplot(VectorToMatrix(XD.init, n, p) ~ Y.logit, main = "Initialisation")
boxplot(VectorToMatrix(XD.hat, n, p) ~ Y.logit, main = "Estimation")


par(mfrow=c(1,1))

plot(log(1 + data$Y[data$Y !=0]), log(1 + pred[data$Y != 0]), main = "Prediction") ; abline(0,1)
plot(log(1 + data$Y[data$Y !=0]), log(1 + pred.init[data$Y != 0]), main = "Prediction") ; abline(0,1)

plot( data$Y[data$Y !=0], pred[data$Y != 0], main = "Prediction") ; abline(0,1)
plot( data$Y[data$Y !=0], pred.init[data$Y != 0], main = "Prediction") ; abline(0,1)



