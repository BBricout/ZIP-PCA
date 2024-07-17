# Test with/without miss

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
# seed <- .Random.seed
library(nloptr)
setwd('/home/robin/RECHERCHE/ECOLOGIE/B-Bricout/ZIP-PCA/codesSR/tests/')
source('../Functions/FunctionsZIP.R'); 
source('../Functions/FunctionsUtils.R'); 
source('../Functions/FunctionsZIPLNmiss.R')
source('../Functions/FunctionsZIPLNmissVec.R')
source('../Functions/FunctionsBB2SR.R')
setwd('/home/robin/RECHERCHE/ECOLOGIE/B-Bricout/ZIP-PCA/')
source('codesBB/FunctionsBB.R')
source('codesBB/UtilsBB.R')
setwd('/home/robin/RECHERCHE/ECOLOGIE/B-Bricout/ZIP-PCA/codesSR/tests/')
simDir <- '../../simulSR/'
resDir <- './'

# Parms
n <- 100; d <- 5; p <- 10; q <- 2; 
seedList <- 1:100; obsList <- c(1, 0.8, 0.5)
# iMiss <- 1; jMiss <- 5; ijMiss <- c(iMiss, jMiss)
load('../../simulSR/ZiPLNsim-n100-d5-p10-q2-seed1-obs80.Rdata')

# Algo parms
iterMax <- 1e3; tolS <- 1e-6; tolXi <- 0
lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))
lbAbs <- c(rep(-Inf, 2*d + q*(p+n-1)), rep(tolS, (n-1) * q))
algoName <- paste0('-tolS', tolS, '-tolXi', tolXi)
optsMMA <- list("algorithm"="NLOPT_LD_MMA", "xtol_rel"=1.0e-8)
optsLBGFS <- list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8)
optsCCSAQ <- list("algorithm"="NLOPT_LD_CCSAQ", "xtol_rel"=1.0e-8)

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
  #                     q=q, method='L-BFGS-B', control=list(fnscale=-1, maxit=iterMax), lower=lb))
  initLogS <-Init2logS(init)
  mmaLogS <- nloptr(x0=c(Mstep2Theta(initLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initLogS$eStep, n=n, d=d, p=p, q=q)),
                    eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=data, q=q, tolXi=tolXi, opts=optsMMA)
  lbfgsLogS <- nloptr(x0=c(Mstep2Theta(initLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initLogS$eStep, n=n, d=d, p=p, q=q)),
                      eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=data, q=q, tolXi=tolXi, opts=optsLBGFS)
  ccsaqLogS <- nloptr(x0=c(Mstep2Theta(initLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initLogS$eStep, n=n, d=d, p=p, q=q)),
                      eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=data, q=q, tolXi=tolXi, opts=optsCCSAQ)
  gdLogS <- optim(par=c(Mstep2Theta(initLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initLogS$eStep, n=n, d=d, p=p, q=q)),
                  fn=ElboLogSVecThetaPsi, gr=ElboLogSGradVecThetaPsi, data=data, tolXi=tolXi, 
                  q=q, method='BFGS', control=list(fnscale=-1, maxit=iterMax))
  # tab <- c(mma=-mma$objective, lbfgs=-lbfgs$objective, gd=gd$value, 
  #          mmaLogS=-mmaLogS$objective, lbfgsLogS=-lbfgsLogS$objective, gdLogS=gdLogS$value)
  # return(list(init=init, mma=mma, lbfgs=lbfgs, gd=gd, mmaLogS=mmaLogS, lbfgsLogS=lbfgsLogS, gdLogS=gdLogS, tab=tab))
  tab <- c(mmaLogS=-mmaLogS$objective, lbfgsLogS=-lbfgsLogS$objective, ccsaqLogS=-ccsaqLogS$objective, gdLogS=gdLogS$value)
  return(list(init=init, mmaLogS=mmaLogS, lbfgsLogS=lbfgsLogS, ccsaqLogS=ccsaqLogS, gdLogS=gdLogS, tab=tab))
}

# Loop
# elboTab <- array(NA, dim=c(length(seedList), length(obsList), 4, 4), 
#                  dimnames=list(seedList, obsList, c('full', 'miss1', 'missRow', 'absRow'), 
#                                c('mmaLogS', 'lbfgsLogS', 'ccsaqLogS', 'gdLogS')))
resFile <- paste0(resDir, 'CompAlgoGD', algoName, '.Rdata')
if(file.exists(resFile)){
  load(resFile)
}else{
  elboTab <- array(NA, dim=c(length(seedList), length(obsList), 4), 
                   dimnames=list(seedList, obsList, c('mmaLogS', 'lbfgsLogS', 'ccsaqLogS', 'gdLogS')))
}
for(seed in seedList){
# for(seed in 8:length(seedList)){
  cat(seed, ':', sep='')
  if(sum(is.na(elboTab[seed, , ])) > 0){
    for(o in 1:length(obsList)){
      obs <- obsList[o]
      cat('', obs)
      # Data
      simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
      simParms <- paste0(simParmsFull, '-obs', 100*obs)
      simName <- paste0('ZiPLNsim', simParms)
      simFile <- paste0(simDir, simName, '.Rdata')
      load(simFile)
      elboTab[seed, o, ] <- FitAll(data=data, q=q)$tab
      # data$Omega[iMiss, jMiss] <- 1
      # # Full data set
      # dataFull <- data; 
      # # 1 missing observation
      # dataMiss1 <- data; dataMiss1$Omega[iMiss, jMiss] <- 0
      # # 1 missing row
      # dataMissRow <- data; dataMissRow$Omega[iMiss, ] <- 0
      # dataAbsRow <- data; 
      # # 1 absent row
      # dataAbsRow$Y <- data$Y[-iMiss, ]
      # dataAbsRow$logFactY <- data$logFactY[-iMiss, ]
      # dataAbsRow$Omega <- data$Omega[-iMiss, ]
      # dataAbsRow$X <- data$X[-which(data$ij[, 1] %in% iMiss), ]
      # dataAbsRow$ij <- cbind(rep(1:(n-1), p), rep(1:p, each=(n-1)))
      # # Fit
      # elboTab[seed, o, 1, ] <- FitAll(data=dataFull, q=q)$tab
      # elboTab[seed, o, 2, ] <- FitAll(data=dataMiss1, q=q)$tab
      # elboTab[seed, o, 3, ] <- FitAll(data=dataMissRow, q=q)$tab
      # elboTab[seed, o, 4, ] <- FitAll(data=dataAbsRow, q=q)$tab
    }
  }
  cat(', ')
  if(seed > 2){
    par(mfrow=c(length(obsList), 1))
    for(o in 1:length(obsList)){
      # boxplot(elboTab[, o, ], ylim=quantile(elboTab[, o, ], prob=c(0.25, 1), na.rm=TRUE))
      matplot(elboTab[, o, ], xlim=c(1, seed), ylim=quantile(elboTab[, o, ], prob=c(0.25, 1), na.rm=TRUE), type='b')
    }
  }
  save(elboTab, file=paste0(resDir, 'CompAlgoGD', algoName, '.Rdata'))
}
