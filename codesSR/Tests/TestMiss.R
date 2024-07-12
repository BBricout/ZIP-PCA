# Test with/without miss

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
# seed <- .Random.seed
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
iterMax <- 1e3; tolS <- 1e-4; tolXi <- 1e-4
lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))
iMiss <- 5; jMiss <- 5; ijMiss <- c(iMiss, jMiss)

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
# S version
gdFile <- paste0(resDir, simName, '-gd.Rdata')
if(!file.exists(gdFile)){
  gdFull <- try(optim(par=c(Mstep2Theta(initFull$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initFull$eStep, n=n, d=d, p=p, q=q)),
                fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataFull, tolXi=tolXi, 
                q=q, method='L-BFGS-B', control=list(fnscale=-1, trace=2, maxit=iterMax), lower=lb))
  if(length(gdFull)==1){gdFull$value <- NA}
  gdMiss1 <- try(optim(par=c(Mstep2Theta(initMiss1$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMiss1$eStep, n=n, d=d, p=p, q=q)),
                fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataMiss1, tolXi=tolXi, 
                q=q, method='L-BFGS-B', control=list(fnscale=-1, trace=2, maxit=iterMax), lower=lb))
  if(length(gdMiss1)==1){gdMiss1$value <- NA}
  gdMissRow <- try(optim(par=c(Mstep2Theta(initMissRow$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMissRow$eStep, n=n, d=d, p=p, q=q)),
                 fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataMissRow, tolXi=tolXi, 
                 q=q, method='L-BFGS-B', control=list(fnscale=-1, trace=2, maxit=iterMax), lower=lb))
  if(length(gdMissRow)==1){gdMissRow$value <- NA}
  gdAbsRow <- try(optim(par=c(Mstep2Theta(initAbsRow$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initAbsRow$eStep, n=n, d=d, p=p, q=q)),
                 fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataAbsRow, tolXi=tolXi, 
                 q=q, method='L-BFGS-B', control=list(fnscale=-1, trace=2, maxit=iterMax), lower=lb))
  if(length(gdAbsRow)==1){gdAbsRow$value <- NA}
  save(gdFull, gdMiss1, gdMissRow, gdAbsRow, file=gdFile)
}else{load(gdFile)}

# logS version
gdLogSFile <- paste0(resDir, simName, '-gdLogS.Rdata')
if(!file.exists(gdLogSFile)){
  iniFullLogS <-Init2Log2(initFull)
  initMiss1LogS <-Init2Log2(initMiss1)
  initMissRowLogS <-Init2Log2(initMissRow)
  initAbsRowLogS <-Init2Log2(initAbsRow)
  ElboVecThetaPsi(c(Mstep2Theta(initFull$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initFull$eStep, n=n, d=d, p=p, q=q)), data=data, q=q)
  ElboLogSVecThetaPsi(c(Mstep2Theta(iniFullLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(iniFullLogS$eStep, n=n, d=d, p=p, q=q)), data=data, q=q)
  gdLogSFull <- optim(par=c(Mstep2Theta(iniFullLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(iniFullLogS$eStep, n=n, d=d, p=p, q=q)),
                      fn=ElboLogSVecThetaPsi, gr=ElboLogSGradVecThetaPsi, data=dataFull, tolXi=tolXi, 
                      q=q, method='BFGS', control=list(fnscale=-1, trace=2, maxit=iterMax))
  gdLogSMiss1 <- optim(par=c(Mstep2Theta(initMiss1$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMiss1$eStep, n=n, d=d, p=p, q=q)),
                       fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataMiss1, tolXi=tolXi, 
                       q=q, method='BFGS', control=list(fnscale=-1, trace=2, maxit=iterMax))
  gdLogSMissRow <- optim(par=c(Mstep2Theta(initMissRow$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initMissRow$eStep, n=n, d=d, p=p, q=q)),
                         fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataMissRow, tolXi=tolXi, 
                         q=q, method='BFGS', control=list(fnscale=-1, trace=2, maxit=iterMax))
  gdLogSAbsRow <- optim(par=c(Mstep2Theta(initAbsRow$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initAbsRow$eStep, n=n, d=d, p=p, q=q)),
                        fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=dataAbsRow, tolXi=tolXi, 
                        q=q, method='BFGS', control=list(fnscale=-1, trace=2, maxit=iterMax))
  save(gdLogSFull, gdLogSMiss1, gdLogSMissRow, gdLogSAbsRow, file=gdLogSFile)
}else{load(gdLogSFile)}

print(rbind(c(gdFull$value, gdMiss1$value, gdMissRow$value, gdAbsRow$value), 
      c(gdLogSFull$value, gdLogSMiss1$value, gdLogSMissRow$value, gdLogSAbsRow$value)))

#################################################################################
# VEM
#################################################################################
# Full dataset
fullFile <- paste0(resDir, simName, '-full.Rdata')
if(!file.exists(fullFile)){
  fitFull <- VemZiPLN(data=dataFull, init=initFull, iterMax=iterMax, tolXi=tolXi, tolS=tolS, plot=TRUE)
  save(fitFull, file=fullFile)
}else{load(fullFile)}
print(c(fitFull$elbo))
plot(fitFull$elboPath, type='b', pch=20, ylim=quantile(fitFull$elboPath, prob=c(0.1, 1)))
abline(h=c(gdLogSFull$value, gdLogSMiss1$value, gdLogSMissRow$value, gdLogSAbsRow$value), col=1:4)

# One missing data
miss1File <- paste0(resDir, simName, '-miss1.Rdata')
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
miss2File <- paste0(resDir, simName, '-missRow.Rdata')
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
miss3File <- paste0(resDir, simName, '-absRow.Rdata')
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
