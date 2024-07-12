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
n <- 100; d <- 5; p <- 10; q <- 2; seed <- 1; obs <- 1
simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
simParms <- paste0(simParmsFull, '-obs', 100*obs)
simName <- paste0('ZiPLNsim', simParms)
simFile <- paste0(simDir, simName, '.Rdata')
load(simFile)

# Algo parms
iterMax <- 1e3; tolS <- 1e-6; tolXi <- 0
lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))
lbAbs <- c(rep(-Inf, 2*d + q*(p+n-1)), rep(tolS, (n-1) * q))
iMiss <- 5; jMiss <- 5; ijMiss <- c(iMiss, jMiss)
algoName <- paste0(simName, 'tolS', tolS, '-tolXi', tolXi)

# Data sets
dataFull <- dataMiss1 <- dataMissRow <- dataAbsRow <- data; 
dataMiss1$Omega[iMiss, jMiss] <- 0
dataMissRow$Omega[iMiss, ] <- 0
dataAbsRow$Y <- data$Y[-iMiss, ]
dataAbsRow$logFactY <- data$logFactY[-iMiss, ]
dataAbsRow$Omega <- data$Omega[-iMiss, ]
dataAbsRow$X <- data$X[-which(data$ij[, 1] %in% iMiss), ]
dataAbsRow$ij <- cbind(rep(1:(n-1), p), rep(1:p, each=(n-1)))

#################################################################################
# Init
#################################################################################
Init2Log2 <- function(init){init$eStep$S <- log(init$eStep$S); return(init)}
initFull <- InitZiPLN(data=dataFull, q=q, tolXi=tolXi)
initMiss1 <- InitZiPLN(data=dataMiss1, q=q, tolXi=tolXi)
initMissRow <- InitZiPLN(data=dataMissRow, q=q, tolXi=tolXi)
initAbsRow <- InitZiPLN(data=dataAbsRow, q=q, tolXi=tolXi)

#################################################################################
# Gradient ascent
#################################################################################
# nlopt
opts <- list("algorithm"="NLOPT_LD_MMA", "xtol_rel"=1.0e-8)
# S version
gdFile <- paste0(resDir, algoName, '-gd.Rdata')
if(!file.exists(gdFile)){
  nlFull <- nloptr(x0=c(Mstep2Theta(initFull$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initFull$eStep, n=n, d=d, p=p, q=q)),
                   eval_f=NegElboVecThetaPsi, eval_grad_f=NegElboGradVecThetaPsi, data=dataFull, q=q, tolXi=tolXi, lb=lb, opts=opts)
  gdFull <- try(optim(par=c(Mstep2Theta(initFull$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initFull$eStep, n=n, d=d, p=p, q=q)),
                fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataFull, tolXi=tolXi, 
                q=q, method='L-BFGS-B', control=list(fnscale=-1, trace=2, maxit=iterMax), lower=lb))
  if(length(gdFull)==1){gdFull$value <- NA}
  nlMiss1 <- nloptr(x0=c(Mstep2Theta(initMiss1$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMiss1$eStep, n=n, d=d, p=p, q=q)),
                    eval_f=NegElboVecThetaPsi, eval_grad_f=NegElboGradVecThetaPsi, data=dataMiss1, q=q, tolXi=tolXi, lb=lb, opts=opts)
  gdMiss1 <- try(optim(par=c(Mstep2Theta(initMiss1$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMiss1$eStep, n=n, d=d, p=p, q=q)),
                fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataMiss1, tolXi=tolXi, 
                q=q, method='L-BFGS-B', control=list(fnscale=-1, trace=2, maxit=iterMax), lower=lb))
  if(length(gdMiss1)==1){gdMiss1$value <- NA}
  nlMissRow <- nloptr(x0=c(Mstep2Theta(initMissRow$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMissRow$eStep, n=n, d=d, p=p, q=q)),
                      eval_f=NegElboVecThetaPsi, eval_grad_f=NegElboGradVecThetaPsi, data=dataMissRow, q=q, tolXi=tolXi, lb=lb, opts=opts)
  gdMissRow <- try(optim(par=c(Mstep2Theta(initMissRow$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMissRow$eStep, n=n, d=d, p=p, q=q)),
                 fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataMissRow, tolXi=tolXi, 
                 q=q, method='L-BFGS-B', control=list(fnscale=-1, trace=2, maxit=iterMax), lower=lb))
  if(length(gdMissRow)==1){gdMissRow$value <- NA}
  nlAbsRow <- nloptr(x0=c(Mstep2Theta(initAbsRow$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initAbsRow$eStep, n=n, d=d, p=p, q=q)),
                     eval_f=NegElboVecThetaPsi, eval_grad_f=NegElboGradVecThetaPsi, data=dataAbsRow, q=q, tolXi=tolXi, lb=lbAbs, opts=opts)
  gdAbsRow <- try(optim(par=c(Mstep2Theta(initAbsRow$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initAbsRow$eStep, n=n, d=d, p=p, q=q)),
                 fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataAbsRow, tolXi=tolXi, 
                 q=q, method='L-BFGS-B', control=list(fnscale=-1, trace=2, maxit=iterMax), lower=lbAbs))
  if(length(gdAbsRow)==1){gdAbsRow$value <- NA}
  save(nlFull, nlMiss1, nlMissRow, nlAbsRow, gdFull, gdMiss1, gdMissRow, gdAbsRow, file=gdFile)
}else{load(gdFile)}

# logS version
gdLogSFile <- paste0(resDir, algoName, '-gdLogS.Rdata')
if(!file.exists(gdLogSFile)){
  initFullLogS <-Init2Log2(initFull); initMiss1LogS <- Init2Log2(initMiss1)
  initMissRowLogS <-Init2Log2(initMissRow); initAbsRowLogS <-Init2Log2(initAbsRow)
  nlLogSFull <- nloptr(x0=c(Mstep2Theta(initFullLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initFullLogS$eStep, n=n, d=d, p=p, q=q)),
                       eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=dataFull, q=q, tolXi=tolXi, opts=opts)
  gdLogSFull <- optim(par=c(Mstep2Theta(initFullLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initFullLogS$eStep, n=n, d=d, p=p, q=q)),
                      fn=ElboLogSVecThetaPsi, gr=ElboLogSGradVecThetaPsi, data=dataFull, tolXi=tolXi, 
                      q=q, method='BFGS', control=list(fnscale=-1, trace=2, maxit=iterMax))
  nlLogSMiss1 <- nloptr(x0=c(Mstep2Theta(initMiss1LogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMiss1LogS$eStep, n=n, d=d, p=p, q=q)),
                        eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=dataMiss1, q=q, tolXi=tolXi, opts=opts)
  gdLogSMiss1 <- optim(par=c(Mstep2Theta(initMiss1$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMiss1$eStep, n=n, d=d, p=p, q=q)),
                       fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataMiss1, tolXi=tolXi, 
                       q=q, method='BFGS', control=list(fnscale=-1, trace=2, maxit=iterMax))
  nlLogSMissRow <- nloptr(x0=c(Mstep2Theta(initMissRowLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMissRowLogS$eStep, n=n, d=d, p=p, q=q)),
                       eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=dataMissRow, q=q, tolXi=tolXi, opts=opts)
  gdLogSMissRow <- optim(par=c(Mstep2Theta(initMissRow$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMissRow$eStep, n=n, d=d, p=p, q=q)),
                         fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataMissRow, tolXi=tolXi, 
                         q=q, method='BFGS', control=list(fnscale=-1, trace=2, maxit=iterMax))
  nlLogSAbsRow <- try(nloptr(x0=c(Mstep2Theta(initAbsRowLogS$mStep, n=n-1, d=d, p=p, q=q), Estep2Psi(initAbsRowLogS$eStep, n=n-1, d=d, p=p, q=q)),
                             eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=dataAbsRow, q=q, tolXi=tolXi, opts=opts))
  if(length(nlLogSAbsRow)==1){nlLogSAbsRow$objective <- NA}
  gdLogSAbsRow <- optim(par=c(Mstep2Theta(initAbsRow$mStep, n=n-1, d=d, p=p, q=q), Estep2Psi(initAbsRow$eStep, n=n-1, d=d, p=p, q=q)),
                        fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataAbsRow, tolXi=tolXi, 
                        q=q, method='BFGS', control=list(fnscale=-1, trace=2, maxit=iterMax))
  save(nlLogSFull, nlLogSMiss1, nlLogSMissRow, nlLogSAbsRow, gdLogSFull, gdLogSMiss1, gdLogSMissRow, gdLogSAbsRow, file=gdLogSFile)
}else{load(gdLogSFile)}

print(rbind(
  c(nlFull$objective, nlMiss1$objective, nlMissRow$objective, nlAbsRow$objective), 
  c(gdFull$value, gdMiss1$value, gdMissRow$value, gdAbsRow$value), 
  c(nlLogSFull$objective, nlLogSMiss1$objective, nlLogSMissRow$objective, nlLogSAbsRow$objective),
  c(gdLogSFull$value, gdLogSMiss1$value, gdLogSMissRow$value, gdLogSAbsRow$value))
)

par(mfrow=c(2, 2))
if(!is.null(gdFull$par)){plot(gdFull$par, gdLogSFull$par); abline(a=0, b=1, h=0, v=0)}
if(!is.null(gdMiss1$par)){plot(gdMiss1$par, gdLogSMiss1$par); abline(a=0, b=1, h=0, v=0)}
if(!is.null(gdMissRow$par)){plot(gdMissRow$par, gdLogSMissRow$par); abline(a=0, b=1, h=0, v=0)}
if(!is.null(gdAbsRow$par)){plot(gdAbsRow$par, gdLogSAbsRow$par); abline(a=0, b=1, h=0, v=0)}

#################################################################################
# VEM
#################################################################################
par(mfrow=c(1, 1))
# Full dataset
fullFile <- paste0(resDir, algoName, '-full.Rdata')
if(!file.exists(fullFile)){
  fitFull <- VemZiPLN(data=dataFull, init=initFull, iterMax=iterMax, tolXi=tolXi, tolS=tolS, plot=TRUE)
  save(fitFull, file=fullFile)
}else{load(fullFile)}
print(c(fitFull$elbo))
plot(fitFull$elboPath, type='b', pch=20, ylim=quantile(fitFull$elboPath, prob=c(0.1, 1)))
abline(h=c(gdLogSFull$value, gdLogSMiss1$value, gdLogSMissRow$value, gdLogSAbsRow$value), col=1:4)

# One missing data
miss1File <- paste0(resDir, algoName, '-miss1.Rdata')
if(!file.exists(miss1File)){
  fitMiss1 <- VemZiPLN(data=dataMiss1, init=initMiss1, iterMax=iterMax, tolXi=tolXi, tolS=tolS, plot=TRUE)
  save(dataMiss1, initMiss1, fitMiss1, file=miss1File)
}else{load(miss1File)}
print(c(fitFull$elbo, fitMiss1$elbo))
plot(fitFull$elboPath, xlim=c(1, max(c(fitFull$iter, fitMiss1$iter))), type='b', pch=20,
     ylim=quantile(c(fitFull$elboPath, fitMiss1$elboPath), prob=c(0.1, 1)))
points(fitMiss1$elboPath, type='b', pch=20, col=2)
abline(h=c(gdLogSFull$value, gdLogSMiss1$value, gdLogSMissRow$value, gdLogSAbsRow$value), col=1:4)

# One missing row
miss2File <- paste0(resDir, algoName, '-missRow.Rdata')
if(!file.exists(miss2File)){
  fitMissRow <- VemZiPLN(data=dataMissRow, init=initMissRow, iterMax=iterMax, tolXi=tolXi, tolS=tolS, plot=TRUE)
  save(dataMissRow, initMissRow, fitMissRow, file=miss2File)
}else{load(miss2File)}
print(c(fitFull$elbo, fitMiss1$elbo, fitMissRow$elbo))
plot(fitFull$elboPath, xlim=c(1, max(c(fitFull$iter, fitMiss1$iter, fitMissRow$iter))), type='b', pch=20,
     ylim=quantile(c(fitFull$elboPath, fitMiss1$elboPath, fitMissRow$elboPath), prob=c(0.1, 1)))
points(fitMiss1$elboPath, type='b', pch=20, col=2)
points(fitMissRow$elboPath, type='b', pch=20, col=3)
abline(h=c(gdLogSFull$value, gdLogSMiss1$value, gdLogSMissRow$value, gdLogSAbsRow$value), col=1:4)

# One absent row
miss3File <- paste0(resDir, algoName, '-absRow.Rdata')
if(!file.exists(miss3File)){
  fitAbsRow <- VemZiPLN(data=dataAbsRow, init=initAbsRow, iterMax=iterMax, tolXi=tolXi, tolS=tolS, plot=TRUE)
  save(dataAbsRow, initAbsRow, fitAbsRow, file=miss3File)
}else{load(miss3File)}
print(c(fitFull$elbo, fitMiss1$elbo, fitMissRow$elbo, fitAbsRow$elbo))
plot(fitFull$elboPath, xlim=c(1, max(c(fitFull$iter, fitMiss1$iter, fitMissRow$iter, fitAbsRow$iter))), type='b', pch=20, 
     ylim=quantile(c(fitFull$elboPath, fitMiss1$elboPath, fitMissRow$elboPath, fitAbsRow$elboPath), prob=c(0.1, 1)))
points(fitMiss1$elboPath, type='b', pch=20, col=2)
points(fitMissRow$elboPath, type='b', pch=20, col=3)
points(fitAbsRow$elboPath, type='b', pch=20, col=4)
abline(h=c(gdFull$value, gdMiss1$value, gdMissRow$value, gdAbsRow$value), col=1:4)
abline(h=c(gdLogSFull$value, gdLogSMiss1$value, gdLogSMissRow$value, gdLogSAbsRow$value), col=1:4, lty=2)

print(rbind(c(gdFull$value, gdMiss1$value, gdMissRow$value, gdAbsRow$value), 
      c(gdLogSFull$value, gdLogSMiss1$value, gdLogSMissRow$value, gdLogSAbsRow$value), 
      c(fitFull$elbo, fitMiss1$elbo, fitMissRow$elbo, fitAbsRow$elbo)))

# ELBO(data=data2, mStep=init2$mStep, eStep=init2$eStep)
# init23 <- init2
# init23$eStep$M <- init2$eStep$M[-iMiss, ]
# init23$eStep$S <- init2$eStep$S[-iMiss, ]
# ELBO(data=data3, mStep=init3$mStep, eStep=init3$eStep)
# ELBO(data=data3, mStep=init23$mStep, eStep=init23$eStep)
