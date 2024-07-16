rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

setwd("~/Documents/ZIP-PCA")

library(Rcpp)
library(PLNmodels)
library(missForest)
library(nloptr)

#-------------------------------------------------------------------------------
# Codes

source("codesBB/FunctionsBB.R")
source("codesBB/UtilsBB.R")

#-------------------------------------------------------------------------------
# Données

load(file='ourfile.Rdata')

Y <- data$Y
X <- data$X

n <- nrow(Y)
p <- ncol(Y)
q <- ncol(true$eStep$M)
d <- ncol(X)

tolXi <- 1e-04 ; tolS <- 1e-04
lBound <- c(rep(-Inf, (2*d)+(p*q)+(n*q)), rep(tolS, n*q))
config <- PLNPCA_param()$config_optim
config$lower_bounds <- lBound

Y.na <- Y
Y.na[1,3] <- NA
data.na <- data
data.na$Y.na <- Y.na
data.na$R <- data.na$Omega <- ifelse(is.na(Y.na), 0, 1)
# data.na$Y <- ifelse(data.na$R == 0, 0, Y)


#-------------------------------------------------------------------------------
# Initialisation des paramètres

init <- Init_ZIP(Y, X, q)
init$M <- matrix(0, n, q)
init.na <- Init_ZIP(data.na$Y.na, X, q)
init.na$M <- matrix(0, n, q)

mStep <- list(gamma = init$D, beta = init$B, C = init$C)
eStep <- list(M = init$M, S = init$S)

eStep$xi <- Elbo_grad_Rcpp(data.na, init, tolXi)$xi


# params <- list(B=as.matrix(mStep$beta),
#                D=as.matrix(mStep$gamma),
#                C=as.matrix(mStep$C),
#                M=as.matrix(eStep$M),
#                S=as.matrix(eStep$S))

parmsInit <- c(mStep$gamma, mStep$beta, as.vector(mStep$C),
               as.vector(t(eStep$M)), as.vector(t(eStep$S)))
parmsLogSInit <- c(mStep$gamma, mStep$beta, as.vector(mStep$C),
                   as.vector(t(eStep$M)), log(as.vector(t(eStep$S))))

parmsInit.na <- c(init.na$D, init.na$B, init$C,
                  as.vector(t(init.na$M)), as.vector(t(init.na$S)))
parmsLogSInit.na <- c(init.na$D, init.na$B, as.vector(init.na$C),
                   as.vector(t(init.na$M)), log(as.vector(t(init.na$S))))

#-------------------------------------------------------------------------------
# Avec ou sans une données

# ElboB <- Elbo(data, init, tolXi)
# ElboB.na <- Elbo_grad_Rcpp(data.na, init, tolXi)$objective
# ElboS.na <- ELBO(data.na, mStep = mStep, eStep = eStep)
# print(c(ElboB, ElboB.na))
# GradB <- Grad(data, init, tolXi)
# GradB.na <- Grad(data.na, init, tolXi)
# 
# # par(mfrow = c(2, 3))
# 
# plot(GradB$gradB, GradB.na$gradB) ; abline(0,1)
# plot(GradB$gradD, GradB.na$gradD) ; abline(0,1)
# plot(GradB$gradC, GradB.na$gradC) ; abline(0,1)
# plot(GradB$gradM, GradB.na$gradM) ; abline(0,1)
# plot(GradB$gradS, GradB.na$gradS) ; abline(0,1)
#-------------------------------------------------------------------------------
# Fonction optimisation

Bobj <- function(parms, data, tolXi){
  params <- Parms2Params(parms=parms, data=data)
  Elbo(data=data, params=params, tolXi)$objective
}
Bobj.neg <- function(parms, data, tolXi){
  -Bobj(parms, data, tolXi)
}
Bgrad <- function(parms, data, tolXi){
  params <- Parms2Params(parms=parms, data=data)
  grads <- Grad(data=data, params=params, tolXi)
  c(as.vector(grads$gradD), as.vector(grads$gradB), as.vector(grads$gradC), 
    as.vector(t(grads$gradM)), as.vector(t(grads$gradS)))
}
Bgrad.neg <- function(parms, data, tolXi){
  -Bgrad(parms, data, tolXi)
}
BobjLogS <- function(parmsLogS, data, tolXi){
  params <- Parms2Params(parms=parmsLogS, data=data)
  params$logS <- params$S
  Elbo_grad_logS_Rcpp(data, params, tolXi)$objective
}
BobjLogS.neg <- function(parmsLogS, data, tolXi){
  -BobjLogS(parmsLogS, data, tolXi)
}
BgradLogS <- function(parmsLogS, data, tolXi){
  params <- Parms2Params(parms=parmsLogS, data=data)
  params$logS <- params$S
  elbo <- Elbo_grad_logS_Rcpp(data, params, tolXi)
  c(as.vector(elbo$gradD), as.vector(elbo$gradB), as.vector(elbo$gradC),
    as.vector(t(elbo$gradM)), as.vector(t(elbo$gradS))*as.vector(t(params$S)))
}
BgradLogS.neg <- function(parmsLogS, data, tolXi){
  -BgradLogS(parmsLogS, data, tolXi)
}


#-------------------------------------------------------------------------------
# Tests avec mon algorithme

# Données complètes
Bfit <- Miss.ZIPPCA(Y, X, q, params = init, config = config)

# Données incomplètes
Bfit.na <- Miss.ZIPPCA(Y.na, X, q, params = init.na, config = config)

# Données complètes avec logS
BfitLogS <- Miss.ZIPPCA.logS(Y, X, params = init)

# Données incomplètes avec logS
BfitLogS.na <- Miss.ZIPPCA(Y.na, X, params = init.na)


#-------------------------------------------------------------------------------
# Tests optimisation

# Données complètes
Optfit <- optim(par=parmsInit, fn=Bobj, gr=Bgrad, data=data,
              method='L-BFGS-B', control=list(fnscale=-1), lower=lBound, tolXi = tolXi)


# Données incomplètes

Optfit.na <- optim(par=parmsInit.na, fn=Bobj, gr=Bgrad, data=data.na,
                   method='L-BFGS-B', control=list(fnscale=-1), lower=lBound, tolXi = tolXi)

# Données complètes avec logS

OptfitLogS <- optim(par=parmsLogSInit, fn=BobjLogS, gr=BgradLogS, data=data,
                  method='BFGS', control=list(fnscale=-1), tolXi = tolXi)

# Données incomplètes avec logS

OptfitLogS.na <- optim(par=parmsLogSInit.na, fn=BobjLogS, gr=BgradLogS, data=data.na,
                  method='BFGS', control=list(fnscale=-1), tolXi = tolXi)

#-------------------------------------------------------------------------------
# Nlopt

opts <- list("algorithm"="NLOPT_LD_CCSAQ",
             "xtol_rel"=1.0e-6)

NLfit <- nloptr(x0 = parmsInit, data = data, eval_f = Bobj.neg, eval_grad_f = Bgrad.neg,
                lb = lBound, tolXi = tolXi, opts = opts)

NLfit.na <- nloptr(x0 = parmsInit.na, data = data.na, eval_f = Bobj.neg, eval_grad_f = Bgrad.neg,
                lb = lBound, tolXi = tolXi, opts = opts)

NLfitLogS <- nloptr(x0 = parmsLogSInit, data = data, eval_f = BobjLogS.neg, eval_grad_f = BgradLogS.neg,
                   tolXi = tolXi, opts = opts)

NLfitLogS.na <- nloptr(x0 = parmsLogSInit.na, data = data.na, eval_f = BobjLogS.neg, eval_grad_f = BgradLogS.neg,
                    tolXi = tolXi, opts = opts)


#-------------------------------------------------------------------------------
# ELBO

# Données complètes
print(c(Bfit$elbo, BfitLogS$elbo, NLfit$objective, NLfitLogS$objective,
        Optfit$value, OptfitLogS$value))

# Données incomplètes
print(c(Bfit.na$elbo, BfitLogS.na$elbo, NLfit.na$objective, NLfitLogS.na$objective,
        Optfit.na$value, OptfitLogS.na$value))
















