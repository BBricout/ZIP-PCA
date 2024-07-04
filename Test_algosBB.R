rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

setwd("~/Documents/ZIP-PCA") # A changer
library(Rcpp)
library(PLNmodels)

# seed <- .Random.seed
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIP.R')
source("codesSR/Functions/FunctionsZIPLNmiss.R")
source("codesSR/Functions/FunctionsZIPLNmissVec.R")
source("codesBB/FunctionsBB.R")
source("codesBB/UtilsBB.R")
source("Data_create.R")


init <- InitZiPLN(data,q) 
params <- list(B = as.matrix(init$mStep$beta), 
               D = as.matrix(init$mStep$gamma), 
               C = as.matrix(init$mStep$C), 
               M = as.matrix(init$eStep$M), 
               S = as.matrix(init$eStep$S))

# plot(init$mStep$beta, params$B) ; abline(0,1)
# plot(init$mStep$gamma, params$D) ; abline(0,1)
# plot(init$mStep$C, params$C) ; abline(0,1)
# plot(init$eStep$M, params$M) ; abline(0,1)
# plot(init$eStep$S, params$S) ; abline(0,1)


# Vérification gradients et ELBO

mStep <- init$mStep ; eStep <- init$eStep

# Belbo_grad <- ElboB(data, params, tolXi)
# 
# Selbo <- ELBO(data=data, mStep=mStep, eStep=eStep)
# SgradS <- matrix(ElboGradVecS(Svec=as.vector(t(eStep$S)), data=data, mStep=mStep, eStep=eStep),n, q, byrow=TRUE)
# SgradM <- matrix(ElboGradVecM(Mvec=as.vector(t(eStep$M)), data=data, mStep=mStep, eStep=eStep),n, q, byrow=TRUE)
# SgradBeta <- ElboGradBeta(beta=mStep$beta, data=data, mStep=mStep, eStep=eStep)
# SgradGamma <- ElboGradGamma(gamma=mStep$gamma, data=data, mStep=mStep, eStep=eStep)
# SgradC <- as.matrix(ElboGradC(vecC=as.vector(mStep$C), data=data, mStep=mStep, eStep=eStep),p,q)
# 
# Belbo_grad$objective ; Selbo
# 
# plot(init$eStep$xi, Belbo_grad$xi) ; abline(0,1)
# 
# plot(Belbo_grad$gradB, SgradBeta) ; abline(0,1)
# plot(Belbo_grad$gradD, SgradGamma) ; abline(0,1)
# plot(Belbo_grad$gradC, SgradC) ; abline(0,1)
# plot(Belbo_grad$gradM, SgradM) ; abline(0,1)
# plot(Belbo_grad$gradS, SgradS) ; abline(0,1)


## Comparaison avec optim


###################################################################################################
## optim_rank_ZIP

# source("codesBB/FunctionsBB.R")


lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(1e-06, n * q))

config <- PLNPCA_param()$config_optim ; tolXi <- 1e-04
config$algorithm <- "MMA" # Par défaut dans config c'est "CCSAQ"
config$lower_bounds <- lb # Si tu veux voir les résultats sans donner la lb il ne faut pas la mettre dans config

out <- Miss.ZIPPCA(Y = data$Y, X = data$X, q, params = params, config = config)

# vem <- VemZiPLN(data=data, init=init, iterMax=5e3)


# Estimations
B.hat <- out$mStep$beta
D.hat <- out$mStep$gamma
C.hat <- out$mStep$C
M.hat <- out$eStep$M
S.hat <- out$eStep$S

XB <- VectorToMatrix(data$X %*% true$mStep$beta, n, p)
XB.hat <- VectorToMatrix(data$X %*% B.hat, n, p)
XD <- VectorToMatrix(data$X %*% true$mStep$gamma, n, p)
XD.hat <- VectorToMatrix(data$X %*% D.hat, n, p)

pred <- exp(XB.hat + M.hat %*% t(C.hat) + 0.5 * (S.hat*S.hat) %*% t(C.hat * C.hat))

Y.logit <- ifelse(data$Y == 0, 0, 1)

# Plot 

par(mfrow=c(2, 3))

plot(out$elboPath, main = "ELBO path")
plot(XD, XD.hat, main = "Estimation de la logistique") ; abline(0,1)
plot(XB, XB.hat, main = "Estimation des régresseurs") ; abline(0,1)
boxplot(VectorToMatrix(XD, n, p) ~ Y.logit, main = "True")
boxplot(VectorToMatrix(XD_ref, n, p) ~ Y.logit, main = "Estimation")
plot(log(1 + data$Y[data$Y !=0]), log(1 + pred[data$Y != 0]), main = "Prediction") ; abline(0,1)


#####################################################################################

# parmsInit <- c(mStep$gamma, mStep$beta, as.vector(mStep$C), 
#                as.vector(t(eStep$M)), as.vector(t(eStep$S)))
# 
# Bobj <- function(parms, data){
#   params <- Parms2Params(parms=parms, data=data)
#   ElboB(data=data, params=params, tolXi = 1e-04)$obj
# }
# Bgrad <- function(parms, data){
#   params <- Parms2Params(parms=parms, data=data)
#   elbo <- ElboB(data=data, params=params, tolXi = 1e-04)
#   c(as.vector(elbo$gradD), as.vector(elbo$gradB), as.vector(elbo$gradC), 
#     as.vector(t(elbo$gradM)), as.vector(t(elbo$gradS)))
# }
# 
# Bfit <- optim(par=parmsInit, fn=Bobj, gr=Bgrad, data=data, 
#               method='L-BFGS-B', control=list(fnscale=-1), lower = lb)
# 
# Bfit$value













# 
# Selbo <- ELBOSophie(data, init$mStep, init$eStep)
# Belbo <- ElboB(data, params, tolXi=tolXi)
# 
# plot(init$eStep$xi, Belbo$xi) ; abline(0,1)
# plot(ComputeXi(data, init$mStep, init$eStep, tolXi=tolXi), Belbo$xi) ; abline(0,1)
# 
# cbind(ELBOSophie(data = data, init$mStep,  init$eStep), unlist(ElboB(data,params, tolXi)))
# Selbo <- ELBO(data, init$mStep, init$eStep)
# 
# 
# ## optim avec lowerbound
# 
# lb1 <- rep(-Inf, 2*d + q*(p+n) + 2)
# lb2 <- rep(1e-04, n * q)
# lb <- c(lb1, lb2)
# sourceCpp("src/optim_rank_ZIP_lb.cpp")
# 
# config$lower_bounds <- lb
# out_lb <- nlopt_optimize_ZIP_lb(data, params, config, tolXi)
# out_lb$objective_values[length(out_lb$objective_values)]
# out_lb$objective_values











