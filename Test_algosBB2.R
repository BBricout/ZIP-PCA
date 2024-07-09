rm(list=ls()); palette('R3')

library(Rcpp)
library(PLNmodels)
library(missForest)
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIP.R')
source("codesSR/Functions/FunctionsZIPLNmiss.R")
source("codesSR/Functions/FunctionsZIPLNmissVec.R")
source("codesBB/FunctionsBB.R")
source("codesBB/UtilsBB.R")

# Data
# seed <- 2; set.seed(seed)
q <- 2
data <- SimZiPLN(n=100, p=10, d=5, q=q)
load(file='ourfile.Rdata')
data$Y.na <- prodNA(data$Y, 0.33)
data$R <- data$Omega <- ifelse(is.na(data$Y.na), 0, 1)
data$Y <- data$Y*data$R
n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X)

# Function
Fit2Parms <- function(mStep, eStep){
  list(B = as.matrix(mStep$beta), D = as.matrix(mStep$gamma), C = as.matrix(mStep$C), 
       M = as.matrix(eStep$M), S = as.matrix(eStep$S))
}

# Init
tolXi <- 1e-06; initS <- 1e-01
init <- InitZiPLN(data, q, tolXi, initS) 
parmInit <- Fit2Parms(mStep=init$mStep, eStep=init$eStep)

# BB
elboB <- ElboB(data=data, params=parmInit, tolXi=tolXi)
# SR
mStep <- init$mStep ; eStep <- init$eStep ;  
Selbo <- ELBO(data=data, mStep=mStep, eStep=eStep)
SgradS <- matrix(ElboGradVecS(Svec=as.vector(t(eStep$S)), data=data, mStep=mStep, eStep=eStep),n, q, byrow=TRUE)
SgradM <- matrix(ElboGradVecM(Mvec=as.vector(t(eStep$M)), data=data, mStep=mStep, eStep=eStep),n, q, byrow=TRUE)
SgradBeta <- ElboGradBeta(beta=mStep$beta, data=data, mStep=mStep, eStep=eStep)
SgradGamma <- ElboGradGamma(gamma=mStep$gamma, data=data, mStep=mStep, eStep=eStep)
SgradC <- as.matrix(ElboGradC(vecC=as.vector(mStep$C), data=data, mStep=mStep, eStep=eStep),p,q)

# Compare
par(mfrow=c(3, 2), pch=20); 
plot(init$eStep$xi, elboB$xi); abline(0,1)
plot(elboB$gradB, SgradBeta); abline(0,1)
plot(elboB$gradD, SgradGamma); abline(0,1)
plot(elboB$gradC, SgradC); abline(0,1)
plot(elboB$gradM, SgradM); abline(0,1)
plot(elboB$gradS, SgradS); abline(0,1)
print(c(elboB$objective, Selbo))
sum(data$R==0)/prod(dim(data$R))

# Optim BB
tolS <- 1e-04
config <- PLNPCA_param()$config_optim
config$algorithm <- "MMA" # Par défaut dans config c'est "CCSAQ"
lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))
config$lower_bounds <- lb # Si tu veux voir les résultats sans donner la lb il ne faut pas la mettre dans config
iterMax <- config$maxeval <- 105
fitBB1 <- Miss.ZIPPCA(Y=data$Y.na, X=data$X, q, params=parmInit, config=config, tolXi=tolXi)
iterMax <- config$maxeval <- 1e3
fitBB2 <- Miss.ZIPPCA(Y=data$Y.na, X=data$X, q, params=parmInit, config=config, tolXi=tolXi)
plot(fitBB1$elboPath, log='y', ylim=c(1, 1e5))

# # Pb divergence
# parmFit1 <- Fit2Parms(mStep=fitBB1$mStep, eStep=fitBB1$eStep)
# parmFit2 <- Fit2Parms(mStep=fitBB2$mStep, eStep=fitBB2$eStep)
# par(mfrow=c(3, 2))
# for(i in 1:length(parmFit1)){plot(parmFit1[[i]], parmFit2[[i]]); abline(0,1)}
# elboB1 <- ElboB(data=data, params=parmFit1, tolXi=tolXi)
# sapply(1:length(elboB1), function(i){sum(is.na(elboB1[[i]]))})
# elboB2 <- ElboB(data=data, params=parmFit2, tolXi=tolXi)
# sapply(1:length(elboB2), function(i){sum(is.na(elboB2[[i]]))})
# print(c(elboB$objective, elboB1$objective, elboB2$objective))

# # Optim SR
# theta <- Mstep2Theta(mStep, n=n, d=d, p=p, q=q)
# psi <- Estep2Psi(eStep, n=n, d=d, p=p, q=q)
# iterMax <- config$maxeval <- 5e2
# fitSR <- optim(par=c(theta, psi), fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=data,
#                q=q, method='L-BFGS-B', control=list(fnscale=-1, trace=2, maxit=iterMax),
#                lower=lb)
# print(c(elboB$objective, elboB1$objective, elboB2$objective, fitSR$value))
# thetaPsi <- fitSR$par
# theta <- thetaPsi[1:((2*d)+(p*q))]
# mStep <- Theta2Mstep(theta, n=n, d=d, p=p, q=q)
# psi <- thetaPsi[(2*d)+(p*q) + (1:(2*(n*q)))]
# eStep <- Psi2Estep(psi, n=n, d=d, p=p, q=q)
# eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=tolXi)
# parmFitSR <- Fit2Parms(mStep=mStep, eStep=eStep)
# par(mfrow=c(3, 2))
# for(i in 1:length(parmFit1)){plot(parmFitSR[[i]], parmFit2[[i]]); abline(0,1)}
