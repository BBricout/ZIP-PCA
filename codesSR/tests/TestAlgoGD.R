# Test with/without miss

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
# seed <- .Random.seed
library(nloptr)
source('../Functions/FunctionsZIP.R'); 
source('../Functions/FunctionsUtils.R'); 
source('../Functions/FunctionsZIPLNmiss.R')
source('../Functions/FunctionsZIPLNmissVec.R')
simDir <- '../../simulSR/'
resDir <- './'

# Data
n <- 100; d <- 5; p <- 10; q <- 2; seed <- 1; obs <- 0.8
simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
simParms <- paste0(simParmsFull, '-obs', 100*obs)
simName <- paste0('ZiPLNsim', simParms)
simFile <- paste0(simDir, simName, '.Rdata')
load(simFile)

# Algo parms
iterMax <- 1e3; tolS <- 1e-6; tolXi <- 0
lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))
lbAbs <- c(rep(-Inf, 2*d + q*(p+n-1)), rep(tolS, (n-1) * q))
iMiss <- 1; jMiss <- 5; ijMiss <- c(iMiss, jMiss)
while(data$Omega[iMiss, jMiss]==0){iMiss <- iMiss+1}
algoName <- paste0(simName, '-tolS', tolS, '-tolXi', tolXi)

# Full data set
dataFull <- data; 
initFull <- InitZiPLN(data=dataFull, q=q)
ELBO(data=dataFull, mStep=initFull$mStep, eStep=initFull$eStep)

# 1 missing observation
dataMiss1 <- data; 
dataMiss1$Omega[iMiss, jMiss] <- 0
initMiss1 <- InitZiPLN(data=dataMiss1, q=q)
ELBO(data=dataMiss1, mStep=initMiss1$mStep, eStep=initMiss1$eStep)

# 1 missing row
dataMissRow <- data; 
dataMissRow$Omega[iMiss, ] <- 0
initMissRow <- InitZiPLN(data=dataMissRow, q=q)
ELBO(data=dataMissRow, mStep=initMissRow$mStep, eStep=initMissRow$eStep)
gradPsiMissRow <- ElboGradVecPsi(Estep2Psi(initMissRow$eStep, n=n, d=d, p=p, q=q), data=dataMissRow, mStep=initMissRow$mStep)
gradPsiMissRow <- Psi2Estep(psi=gradPsiMissRow, n=n, d=d, p=p, q=q)
head(gradPsiMissRow$M)

# 1 absent row
dataAbsRow <- data; 
dataAbsRow$Y <- data$Y[-iMiss, ]
dataAbsRow$logFactY <- data$logFactY[-iMiss, ]
dataAbsRow$Omega <- data$Omega[-iMiss, ]
dataAbsRow$X <- data$X[-which(data$ij[, 1] %in% iMiss), ]
dataAbsRow$ij <- cbind(rep(1:(n-1), p), rep(1:p, each=(n-1)))
initAbsRow <- initMissRow; 
initAbsRow$eStep$M <- initAbsRow$eStep$M[-iMiss, ]; initAbsRow$eStep$S <- initAbsRow$eStep$S[-iMiss, ]
initAbsRow$eStep$xi <- initAbsRow$eStep$xi[-iMiss, ]
ELBO(data=dataMissRow, mStep=initMissRow$mStep, eStep=initMissRow$eStep)
gradPsiAbsRow <- ElboGradVecPsi(Estep2Psi(initMissRow$eStep, n=n-1, d=d, p=p, q=q), data=dataAbsRow, mStep=initAbsRow$mStep)
gradPsiAbsRow <- Psi2Estep(psi=gradPsiAbsRow, n=n-1, d=d, p=p, q=q)
head(gradPsiAbsRow$M)

# ELBO
par(mfrow=c(2, 2))
plot(ElboGradVecTheta(Mstep2Theta(initMissRow$mStep, n=n, d=d, p=p), data=dataMissRow, eStep=initMissRow$eStep), 
     ElboGradVecTheta(Mstep2Theta(initMissRow$mStep, n=n, d=d, p=p), data=dataAbsRow, eStep=initAbsRow$eStep)); abline(0, 1, h=0, v=1)
plot(gradPsiMissRow$M[-iMiss, ], gradPsiAbsRow$M); abline(0, 1, h=0, v=1)
plot(gradPsiMissRow$S[-iMiss, ], gradPsiAbsRow$S); abline(0, 1, h=0, v=1)

# Functions
Init2logS <- function(init){init$eStep$S <- log(init$eStep$S); return(init)}
FitAll <- function(data, q){
  n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X); 
  init <- InitZiPLN(data=data, q=q, tolXi=tolXi)
  # lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))
  # mma <- try(nloptr(x0=c(Mstep2Theta(init$mStep, n=n, d=d, p=p, q=q), Estep2Psi(init$eStep, n=n, d=d, p=p, q=q)),
  #                   eval_f=NegElboVecThetaPsi, eval_grad_f=NegElboGradVecThetaPsi, data=data, q=q, tolXi=tolXi, lb=lb, opts=optsLBGFS))
  # lbfgs <- try(nloptr(x0=c(Mstep2Theta(init$mStep, n=n, d=d, p=p, q=q), Estep2Psi(init$eStep, n=n, d=d, p=p, q=q)),
  #               eval_f=NegElboVecThetaPsi, eval_grad_f=NegElboGradVecThetaPsi, data=data, q=q, tolXi=tolXi, lb=lb, opts=optsMMA))
  # gd <- try(optim(par=c(Mstep2Theta(init$mStep, n=n, d=d, p=p, q=q), Estep2Psi(init$eStep, n=n, d=d, p=p, q=q)),
  #                     fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=data, tolXi=tolXi, 
  #                     q=q, method='L-BFGS-B', control=list(fnscale=-1, trace=2, maxit=iterMax), lower=lb))
  initLogS <-Init2logS(init)
  mmaLogS <- nloptr(x0=c(Mstep2Theta(initLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initLogS$eStep, n=n, d=d, p=p, q=q)),
                        eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=data, q=q, tolXi=tolXi, opts=optsMMA)
  lbfgsLogS <- nloptr(x0=c(Mstep2Theta(initLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initLogS$eStep, n=n, d=d, p=p, q=q)),
                    eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=data, q=q, tolXi=tolXi, opts=optsLBGFS)
  ccsaqLogS <- nloptr(x0=c(Mstep2Theta(initLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initLogS$eStep, n=n, d=d, p=p, q=q)),
                     eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=data, q=q, tolXi=tolXi, opts=optsCCSAQ)
  gdLogS <- optim(par=c(Mstep2Theta(initLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initLogS$eStep, n=n, d=d, p=p, q=q)),
                      fn=ElboLogSVecThetaPsi, gr=ElboLogSGradVecThetaPsi, data=data, tolXi=tolXi, 
                      q=q, method='BFGS', control=list(fnscale=-1, trace=2, maxit=iterMax))
  # tab <- c(mma=-mma$objective, lbfgs=-lbfgs$objective, gd=gd$value, 
  #          mmaLogS=-mmaLogS$objective, lbfgsLogS=-lbfgsLogS$objective, gdLogS=gdLogS$value)
  # return(list(init=init, mma=mma, lbfgs=lbfgs, gd=gd, mmaLogS=mmaLogS, lbfgsLogS=lbfgsLogS, gdLogS=gdLogS, tab=tab))
  tab <- c(mmaLogS=-mmaLogS$objective, lbfgsLogS=-lbfgsLogS$objective, ccsaqLogS=ccsaqLogS$objective, gdLogS=gdLogS$value)
  return(list(init=init, mmaLogS=mmaLogS, lbfgsLogS=lbfgsLogS, ccsaqLogS=ccsaqLogS, gdLogS=gdLogS, tab=tab))
}

# Fit
optsMMA <- list("algorithm"="NLOPT_LD_MMA", "xtol_rel"=1.0e-8)
optsLBGFS <- list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8)
optsCCSAQ <- list("algorithm"="NLOPT_LD_CCSAQ", "xtol_rel"=1.0e-8)
full <- FitAll(data=dataFull, q=q)
miss1 <- FitAll(data=dataMiss1, q=q)
missRow <- FitAll(data=dataMissRow, q=q)
absRow <- FitAll(data=dataAbsRow, q=q)

res <- rbind(full$tab, miss1$tab, missRow$tab, absRow$tab)
rownames(res) <- c('full', 'miss1', 'missRow', 'absRow')
res
