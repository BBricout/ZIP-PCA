rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

setwd("~/Documents/ZIP-PCA")

library(Rcpp)
library(PLNmodels)
library(missForest)



#--- Codes SR
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIP.R')
source("codesSR/Functions/FunctionsZIPLNmiss.R")
source("codesSR/Functions/FunctionsZIPLNmissVec.R")

#--- Codes BB

# sourceCpp("src/optim_rank_ZIP_new.cpp")
source("codesBB/FunctionsBB.R")
source("codesBB/UtilsBB.R")

#--- Données

load(file='ourfile.Rdata')

n <- nrow(data$Y)
p <- ncol(data$Y)
q <- ncol(true$eStep$M)
d <- ncol(data$X)

tolXi <- 1e-04

#--- Initialisation

init <- InitZiPLN(data,q,tolXi) 
params <- list(B = as.matrix(init$mStep$beta), 
               D = as.matrix(init$mStep$gamma), 
               C = as.matrix(init$mStep$C), 
               M = as.matrix(init$eStep$M), 
               S = as.matrix(init$eStep$S))

#--- Tests sur la fonction Elbo_grad

mStep <- init$mStep ; eStep <- init$eStep

Belbo <- Elbo_grad_Rcpp(data, params, tolXi)
# Selbo <- ELBO(data=data, mStep=init$mStep, eStep=init$eStep)
# Belbo$objective
# Selbo[6]
# 
# 
# SgradS <- matrix(ElboGradVecS(Svec=as.vector(t(eStep$S)), data=data, mStep=mStep, eStep=eStep),n, q, byrow=TRUE)
# SgradM <- matrix(ElboGradVecM(Mvec=as.vector(t(eStep$M)), data=data, mStep=mStep, eStep=eStep),n, q, byrow=TRUE)
# SgradBeta <- ElboGradBeta(beta=mStep$beta, data=data, mStep=mStep, eStep=eStep)
# SgradGamma <- ElboGradGamma(gamma=mStep$gamma, data=data, mStep=mStep, eStep=eStep)
# SgradC <- as.matrix(ElboGradC(vecC=as.vector(mStep$C), data=data, mStep=mStep, eStep=eStep),p,q)
# 
# plot(init$eStep$xi, Belbo$xi) ; abline(0,1)
# plot(Belbo$gradB, SgradBeta) ; abline(0,1)
# plot(Belbo$gradD, SgradGamma) ; abline(0,1)
# plot(Belbo$gradC, SgradC) ; abline(0,1)
# plot(Belbo$gradM, SgradM) ; abline(0,1)
# plot(Belbo$gradS, SgradS) ; abline(0,1)

#--- Configuration

tolS <- 1e-04


lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))

config <- PLNPCA_param()$config_optim
# config$algorithm <- "MMA" # Par défaut dans config c'est "CCSAQ"
config$lower_bounds <- lb

data$Y.na <- prodNA(data$Y, noNA = 0.33)


out <- Miss.ZIPPCA(Y = data$Y.na, X = data$X, q, params = params, config = config)


B.hat <- out$mStep$beta
D.hat <- out$mStep$gamma
C.hat <- out$mStep$C
M.hat <- out$eStep$M
S.hat <- out$eStep$S

XB <- VectorToMatrix(data$X %*% true$mStep$beta, n, p)
XB.hat <- VectorToMatrix(data$X %*% B.hat, n, p)
XB.init <- VectorToMatrix(data$X %*% params$B, n, p)


XD <- VectorToMatrix(data$X %*% true$mStep$gamma, n, p)
XD.hat <- VectorToMatrix(data$X %*% D.hat, n, p)
XD.init <- VectorToMatrix(data$X %*% params$D, n, p)


pred.init <- pred <- exp(XB.init) #  + params$M %*% t(params$C) + 0.5 * (params$S*params$S) %*% t(params$C * params$C))
pred <- exp(XB.hat + M.hat %*% t(C.hat) + 0.5 * (S.hat*S.hat) %*% t(C.hat * C.hat))

Y.logit <- ifelse(data$Y == 0, 0, 1)

#----------------------------  Plot 

#--------------- ELBO 
par(mfrow=c(1, 1))
plot(out$elboPath)
plot(out$elboPath[out$elboPath > 1*Belbo$objective],type='l', main = "ELBO path",ylab='ELBO',xlab='iter')
plot(D.hat, params$D) ; abline(0,1)

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


par(mfrow=c(2,2))

plot(log(1 + data$Y[data$Y !=0]), log(1 + pred[data$Y != 0]), main = "Prediction") ; abline(0,1)
plot(log(1 + data$Y[data$Y !=0]), log(1 + pred.init[data$Y != 0]), main = "Prediction") ; abline(0,1)

plot( data$Y[data$Y !=0], pred[data$Y != 0], main = "Prediction") ; abline(0,1)
plot( data$Y[data$Y !=0], pred.init[data$Y != 0], main = "Prediction") ; abline(0,1)



