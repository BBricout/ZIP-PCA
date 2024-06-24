# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
library(Rcpp)
library(PLNmodels)

# seed <- .Random.seed
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIPLN.R')

source("FunctionsBB.R")

simDir <- 'SimulationsBB/datasim/'
# Parms: many small sims
n <- 100; d <- 5; p <- 10; q <- 2
baseSimName <- 'ZiPLNsim'
seedList <- 1:10; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)

# # Parms: one big sim
# n <- 500; d <- 20; p <- 30; q <- 5
# seedList <- 1; seedNb <- length(seedList)
# obsList <- c(0.6); obsNb <- length(obsList)

#X0 <- NULL
 # Same X for all
seedX <- 1; set.seed(seedX)
#X0 <- matrix(rnorm(n*p*d), n*p, d); X0[, 1] <- 1; 
baseSimName <- paste0(baseSimName, '-sameX', seedX)

# Simul
for(seed in 1:10){
  set.seed(seed)
  simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
  simNameFull <- paste0(baseSimName, simParmsFull)
  simFileFull <- paste0(simDir, simNameFull, '-noMiss.Rdata')
  if(!file.exists(simFileFull)){
    sim <- SimZiPLN(n=n, p=p, d=d, q=q, X=X0)
    obsTresh <- matrix(runif(n*p), n, p)
    data <- list(X=sim$X, Y=sim$Y, obsTresh=obsTresh, ij=sim$ij, logFactY=lgamma(sim$Y+1))
    true <- list(mStep=list(gamma=sim$gamma, beta=sim$beta, C=sim$C),
                 eStep=list(M=sim$W, S=matrix(1e-4, n, q), xi=matrix(plogis(sim$X%*%sim$gamma), n, p)),
                 latent=list(U=sim$U, W=sim$W, Z=sim$Z, Yall=sim$Yall))
    save(data, true, file=simFileFull)
    for(oo in 1:obsNb){
      obs <- obsList[oo]
      simParms <- paste0(simParmsFull, '-obs', 100*obs)
      simName <- paste0(baseSimName, simParms)
      simFile <- paste0(simDir, simName, '.Rdata')
      Omega <- 1*(obsTresh <= obs)
      data <- list(X=sim$X, Y=sim$Y, Omega=Omega, ij=sim$ij, logFactY=lgamma(sim$Y+1))
      save(data, file=simFile)
    }
  }
}




# Fit de Stephane
start_time <- Sys.time()
for(seed in 1:1){
    simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
    simNameFull <- paste0(baseSimName, simParmsFull)
    simFileFull <- paste0(simDir, simNameFull, '-noMiss.Rdata')
    load(file=simFileFull)
    
    
    ################ ajustement SR
    init <- InitZiPLN(data) #! Fichiers différnets pour miss et pas miss
    # vem <- VemZiPLN(data=data, init=init, iterMax=5e3)
    # save(init, vem, file=fitFile)
    ################ ajustement BB
    
    for(oo in 1:obsNb){
        obs <- obsList[oo]
        simParms <- paste0(simParmsFull, '-obs', 100*obs)
        simName <- paste0(baseSimName, simParms)
        simFile <- paste0(simDir, simName, '.Rdata')
      ################ ajustement SR
       
        
    }
    }

end_time <- Sys.time()
execution_time_real <- end_time - start_time

#############################################################
# Fit de Barbara

# start_timeB <- Sys.time()
# outB <- Miss.ZIPPCA(data$Y, data$X, q, params=params)
# end_timeB <- Sys.time()
# execution_time_realB <- end_timeB - start_timeB
# 
# outB$elbo
# vem$elbo
# 
# plot(log(1 + data$Y[data$Y != 0]), log(1 + vem$pred$A[data$Y != 0])) ; abline(0,1)
# plot(log(1 + data$Y[data$Y != 0]), log(1 + outB$pred$A[data$Y != 0])) ; abline(0,1)
# plot(log(1 + vem$pred$A[data$Y != 0]), log(1 + outB$pred$A[data$Y != 0])) ; abline(0,1)
# 
# boxplot(vem$pred$nu[data$Y == 0], vem$pred$nu[data$Y != 0])
# 
# boxplot(outB$pred$nu[data$Y == 0], outB$pred$nu[data$Y != 0])

set.seed(NULL)
data$Omega <- matrix(1, n, p)
data$Omega <- matrix(rbinom(n*p, 1, 0.9), n, p)
data$Y[which(data$Omega==0)] <- 23
data$R <- data$Omega
mStep <- init$mStep; eStep <- init$eStep
mStep$gamma <- rnorm(d); mStep$beta <- rnorm(d); mStep$C <- matrix(rnorm(p*q), p, q)
eStep$M <- matrix(rnorm(n*q), n, q); eStep$S <- matrix(exp(0.01*rnorm(n*q)), n, q);
params <- list(B=as.matrix(mStep$beta), 
               D=as.matrix(mStep$gamma), 
               C=as.matrix(mStep$C), 
               M=as.matrix(eStep$M), 
               S=as.matrix(eStep$S))
eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=0)

source("codesSR/Functions/FunctionsZIPLNmiss.R")
source("codesSR/Functions/FunctionsZIPLNmissVec.R")

Selbo <- ELBO(data=data, mStep=mStep, eStep=eStep)
SgradS <- matrix(ElboGradVecS(Svec=as.vector(t(eStep$S)), data=data, mStep=mStep, eStep=eStep),n, q, byrow=TRUE)
SgradM <- matrix(ElboGradVecM(Mvec=as.vector(t(eStep$M)), data=data, mStep=mStep, eStep=eStep),n, q, byrow=TRUE)
SgradBeta <- ElboGradBeta(beta=mStep$beta, data=data, mStep=mStep, eStep=eStep)
SgradGamma <- ElboGradGamma(gamma=mStep$gamma, data=data, mStep=mStep, eStep=eStep)
SgradC <- as.matrix(ElboGradC(vecC=as.vector(mStep$C), data=data, mStep=mStep, eStep=eStep),p,q)

ElboBfun <- ElboB(data=data, params=params)
ElboBfun$objective
Belbo <- ElboBfun$objective
BgradS <- ElboBfun$gradS
BgradM <-ElboBfun$gradM
BgradBeta <- ElboBfun$gradB
BgradGamma <- ElboBfun$gradD
BgradC <- ElboBfun$gradC

par(mfrow=c(4, 4))
plot(mStep$beta, params$B); abline(0, 1, v=0, h=0)
plot(mStep$gamma, params$D); abline(0, 1, v=0, h=0)
plot(mStep$C, params$C); abline(0, 1, v=0, h=0)
plot(eStep$M, params$M); abline(0, 1, v=0, h=0)
plot(eStep$S, params$S); abline(0, 1, v=0, h=0)
plot(as.vector(data$Y), ElboBfun$vecY); abline(0, 1, v=0, h=0)
plot(as.vector(eStep$xi), ElboBfun$vecxi); abline(0, 1, v=0, h=0)
plot(NuMuA(data=data, mStep=mStep, eStep=eStep)$A, ElboBfun$A); abline(0, 1, v=0, h=0)
plot(NuMuA(data=data, mStep=mStep, eStep=eStep)$nu, ElboBfun$nu); abline(0, 1, v=0, h=0)

plot(BgradBeta, SgradBeta); abline(0, 1, v=0, h=0)
plot(BgradGamma, SgradGamma); abline(0, 1, v=0, h=0)
plot(BgradC, SgradC); abline(0, 1, v=0, h=0)
plot(BgradM, SgradM); abline(0, 1, v=0, h=0)
plot(BgradS, SgradS); abline(0, 1, v=0, h=0)

print(c(Selbo, ElboBfun$objective))




###############################################################################
# Optim
Parms2Steps <- function(parms, data){
  n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X)
  mStep <- list(gamma=parms[1:d], beta=parms[d+(1:d)], 
                C=matrix(parms[(2*d)+(1:(p*q))], p, q))
  eStep <- list(M=matrix(parms[(2*d)+(p*q)+(1:(n*q))], n, q, byrow=TRUE), 
                S=matrix(parms[(2*d)+(p*q)+(n*q)+(1:(n*q))], n, q, byrow=TRUE))
  eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=0)
  return(list(m=mStep, e=eStep))  
}
Sobj <- function(parms, data){
  steps <- Parms2Steps(parms=parms, data=data)
  ELBO(data=data, mStep=steps$m, eStep=steps$e)
}
Sgrad <- function(parms, data){
  steps <- Parms2Steps(parms=parms, data=data)
  c(ElboGradGamma(gamma=steps$m$gamma, data=data, mStep=steps$m, eStep=steps$e),
    ElboGradBeta(beta=steps$m$beta, data=data, mStep=steps$m, eStep=steps$e), 
    ElboGradC(vecC=as.vector(steps$m$C), data=data, mStep=steps$m, eStep=steps$e), 
    ElboGradVecM(Mvec=as.vector(t(steps$e$M)), data=data, mStep=steps$m, eStep=steps$e), 
    ElboGradVecS(Svec=as.vector(t(steps$e$S)), data=data, mStep=steps$m, eStep=steps$e))
}
Parms2Params <- function(parms, data){
  n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X)
  steps <- Parms2Steps(parms=parms, data=data)
  return(list(B=as.matrix(steps$m$beta), 
                 D=as.matrix(steps$m$gamma), 
                 C=as.matrix(steps$m$C), 
                 M=as.matrix(steps$e$M), 
                 S=as.matrix(steps$e$S)))
}
Bobj <- function(parms, data){
  params <- Parms2Params(parms=parms, data=data)
  ElboB(data=data, params=params)$obj
}
Bgrad <- function(parms, data){
  params <- Parms2Params(parms=parms, data=data)
  elbo <- ElboB(data=data, params=params)
  c(as.vector(elbo$gradD), as.vector(elbo$gradB), as.vector(elbo$gradC), 
    as.vector(t(elbo$gradM)), as.vector(t(elbo$gradS)))
}
BobjLogS <- function(parmsLogS, data){
  params <- Parms2Params(parms=parmsLogS, data=data)
  params$S <- exp(params$S)
  ElboB(data=data, params=params)$obj
}
BgradLogS <- function(parmsLogS, data){
  params <- Parms2Params(parms=parmsLogS, data=data)
  params$S <- exp(params$S)
  elbo <- ElboB(data=data, params=params)
  c(as.vector(elbo$gradD), as.vector(elbo$gradB), as.vector(elbo$gradC),
    as.vector(t(elbo$gradM)), as.vector(t(elbo$gradS))*as.vector(t(params$S)))
}

parmsInit <- c(mStep$gamma, mStep$beta, as.vector(mStep$C), 
               as.vector(t(eStep$M)), as.vector(t(eStep$S)))
parmsLogSInit <- c(mStep$gamma, mStep$beta, as.vector(mStep$C), 
                   as.vector(t(eStep$M)), log(as.vector(t(eStep$S))))

tolS <- 1e-12
lBound <- c(rep(-Inf, (2*d)+(p*q)+(n*q)), rep(tolS, n*q))
print(c(Sobj(parms=parmsInit, data=data), Bobj(parms=parmsInit, data=data)))
plot(Sgrad(parms=parmsInit, data=data), Bgrad(parms=parmsInit, data=data)); abline(0, 1, v=0, h=0)
max(abs(Sgrad(parms=parmsInit, data=data) - Bgrad(parms=parmsInit, data=data)))

Sfit <- optim(par=parmsInit, fn=Sobj, gr=Sgrad, data=data, 
              method='L-BFGS-B', control=list(fnscale=-1), lower=lBound)
Sfit$value

# vals <- NULL; i <- 0; path <- matrix(NA, 1e3, length(lBound))
# trace(Bgrad, exit = quote(print(c(returnValue()))))
# trace(Bobj, exit = quote(print(c(i, returnValue(), parms))))
Bfit <- optim(par=parmsInit, fn=Bobj, gr=Bgrad, data=data, 
              method='L-BFGS-B', control=list(fnscale=-1), lower=lBound)
# untrace(Bgrad)
# untrace(Bobj)
Bfit$value

BobjLogS(parmsLogS=parmsLogSInit, data=data)
plot(Bgrad(parms=parmsInit, data=data), 
     BgradLogS(parms=parmsLogSInit, data=data)); abline(0, 1, v=0, h=0)
BfitLogS <- optim(par=parmsLogSInit, fn=BobjLogS, gr=BgradLogS, data=data, 
              method='BFGS', control=list(fnscale=-1))
# untrace(Bgrad)
# untrace(Bobj)
BfitLogS$value

# negBobj <- function(parms, data){-Bobj(parms, data)}
# negBgrad <- function(parms, data){-Bgrad(parms, data)}
# library(nloptr)
# opts <- list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8)
# # opts <- list("algorithm"="NLOPT_LD_CCSAQ", "xtol_rel"=1.0e-8)
# Bfit <- nloptr(x0=parmsInit, eval_f=negBobj, eval_grad_f=negBgrad, data=data, 
#                lb=lBound, opts=opts)
# -Bfit$obj

plot(Sfit$par, Bfit$par); abline(0, 1, v=0, h=0)

sel <- (1:((2*d)+(p*q)+(n*q)))
plot(as.data.frame(cbind(Sfit$par[sel], Bfit$par[sel], BfitLogS$par[sel])))
print(rbind(
  c(Sobj(parms=parmsInit, data=data), Bobj(parms=parmsInit, data=data), 
    BobjLogS(parmsLogS=parmsLogSInit, data=data)),
  c(Sfit$value, Bfit$value, BfitLogS$value)
  ))
