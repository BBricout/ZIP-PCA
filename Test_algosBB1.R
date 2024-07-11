rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

#setwd("~/Documents/ZIP-PCA") # A changer


########################## load data ################
source("Data_create.R")
save(data,true,file="ourfile.Rdata")
rm(list=ls())

load(file='ourfile.Rdata')
################################
# seed <- .Random.seed
library(Rcpp)
library(PLNmodels)
library(missForest)

#--- Codes step
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIP.R')
source("codesSR/Functions/FunctionsZIPLNmiss.R")
source("codesSR/Functions/FunctionsZIPLNmissVec.R")
#----- Codes barbara
source("codesBB/FunctionsBB.R")
source("codesBB/UtilsBB.R")


data$Y.na <- prodNA(data$Y, 0.1)
data$R <- data$Omega <- ifelse(is.na(data$Y.na), 0,1)
# data$Y <- data$Y * data$R


#names(data)
# databis <- list()
# databis$Y <- data$Y[-51,]
# databis$logFactY <- data$logFactY[-51,]
# 
# pastis <- which(data$ij[,1]==51)
# databis$X <- data$X[-pastis,]
# databis$Omega <- data$Omega[-51,]
# databis$R <- databis$Omega
# data<- databis
n <- nrow(data$Y)
p <- ncol(data$Y)
q <- ncol(true$eStep$M)
d <- ncol(data$X)

#data$ij <- matrix(0,n*p,2)
#data$ij[,1] = rep((1:n),p)
#data$ij[,2] = rep((1:p),n)

                      
                      
################### INIT with méthode de stéphane
tolXi <- 1e-04
initS <- 1e-01
#-----------init Stephane
init <- InitZiPLN(data,q,tolXi,initS) 
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


#############################################
# Vérification gradients et ELBO


#---------- Barbara



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


Belbo_grad <- ElboB(data, params, tolXi)

#---------- Stephane
mStep <- init$mStep ; eStep <- init$eStep ;  
Selbo <- ELBO(data=data, mStep=mStep, eStep=eStep)
SgradS <- matrix(ElboGradVecS(Svec=as.vector(t(eStep$S)), data=data, mStep=mStep, eStep=eStep),n, q, byrow=TRUE)
SgradM <- matrix(ElboGradVecM(Mvec=as.vector(t(eStep$M)), data=data, mStep=mStep, eStep=eStep),n, q, byrow=TRUE)
SgradBeta <- ElboGradBeta(beta=mStep$beta, data=data, mStep=mStep, eStep=eStep)
SgradGamma <- ElboGradGamma(gamma=mStep$gamma, data=data, mStep=mStep, eStep=eStep)
SgradC <- as.matrix(ElboGradC(vecC=as.vector(mStep$C), data=data, mStep=mStep, eStep=eStep),p,q)

plot(Belbo_grad$xi, init$eStep$xi)

# unlist(Belbo_grad[-1])
#Selbo

res_ELBO <- c(Belbo_grad$objective , Selbo[6])
names(res_ELBO) = c('B','S')
print(res_ELBO)

plot(init$eStep$xi, Belbo_grad$xi) ; abline(0,1)
plot(Belbo_grad$gradB, SgradBeta) ; abline(0,1)
plot(Belbo_grad$gradD, SgradGamma) ; abline(0,1)
plot(Belbo_grad$gradC, SgradC) ; abline(0,1)
plot(Belbo_grad$gradM, SgradM) ; abline(0,1)
plot(Belbo_grad$gradS, SgradS) ; abline(0,1)



## Comparaison avec optim


###################################################################################################
## optim_rank_ZIP

# source("codesBB/FunctionsBB.R")



lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(1e-06, n * q))

config <- PLNPCA_param()$config_optim
config$algorithm <- "MMA" # Par défaut dans config c'est "CCSAQ"
config$lower_bounds <- lb # Si tu veux voir les résultats sans donner la lb il ne faut pas la mettre dans config
# config$maxeval <- 

out <- Miss.ZIPPCA(Y = data$Y.na, X = data$X, q, params = params, config = config, tolXi)

#vem <- VemZiPLN(data=data, init=init, iterMax=5e3)


# Estimations
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
plot(out$elboPath[out$elboPath< -1*Belbo_grad$objective],type='l', main = "ELBO path",ylab='ELBO',xlab='iter')
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













