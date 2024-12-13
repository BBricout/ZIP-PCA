# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

# seed <- .Random.seed
source('Functions/FunctionsUtils.R')
source('Functions/FunctionsZIP.R')
source('Functions/FunctionsZIPLNmiss.R')
source('Functions/FunctionsZIPLNmissVec.R')
simDir <- '../simulSR/'

# Parms: 
# n <- 100; d <- 5; p <- 10; q <- 2; obs <- 0.70; seed <- 1
n <- 100; d <- 2; p <- 5; q <- 2; coefC <- 0.1; seed <- 3; obs <- 1
simParms <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-coefC', coefC, '-seed', seed, '-obs', 100*obs)
baseSimName <- 'ZiPLNsim'; baseFitName <- 'ZiPLNfit'; 
seedList <- 1:100; seedNb <- length(seedList)
obsList <- c(1, 0.9); obsNb <- length(obsList)

# 1 example
# load(paste0(simDir, simName, simParms, '.Rdata'))
# load(paste0(simDir, fitName, simParms, '.Rdata'))
# vem$eStep$xi <- ComputeXi(data=data, mStep=vem$mStep, eStep=vem$eStep)
# vem$elbo <- ELBO(data=data, mStep=vem$mStep, eStep=vem$eStep)
# vemVec <- VemZiPLNvec(data, init)

# mStep <- init$mStep; theta <- Mstep2Theta(mStep, n=n, d=d, p=p, q=q)
# eStep <- init$eStep; psi <- Estep2Psi(eStep, n=n, d=d, p=p, q=q)
# tolS <- 1e-4; iterMax <- 1e2; tol <- 1e-4; tolXi <- 1e-4
# elboPath <- rep(NA, iterMax)
# iter <- 1; diff <- 2*tol
# elboPath[iter] <- ElboVecThetaPsi(thetaPsi=c(theta, psi), data=data, q=q)
# while((iter & iterMax) & (diff > tol)){
#   iter <- iter+1
#   # # VE + M
#   # fit <- optim(par=c(theta, psi), fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=data,
#   #              q=q, method='L-BFGS-B', control=list(fnscale=-1),
#   #              lower=c(rep(-Inf, ((2*d) + (p*q) + (n*q))), rep(tolS, (n*q))))
#   # thetaNew <- fit$par[1:((2*d)+(p*q))]; psiNew <- fit$par[-(1:((2*d)+(p*q)))]
#   # VE
#   fitVE <- optim(par=psi, fn=ElboVecPsi, gr=ElboGradVecPsi, data=data, mStep=mStep,
#                method='L-BFGS-B', control=list(fnscale=-1),
#                lower=c(rep(-Inf, (n*q)), rep(tolS, (n*q))))
#   psiNew <- fitVE$par; eStepNew <- Psi2Estep(psiNew, n=n, d=d, p=p, q=q)
#   eStepNew$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStepNew, tolXi=tolXi)
#   # M
#   fitM <- optim(par=theta, fn=ElboVecTheta, gr=ElboGradVecTheta, data=data, eStep=eStepNew, 
#                method='BFGS', control=list(fnscale=-1))
#   thetaNew <- fitM$par; mStepNew <- Theta2Mstep(thetaNew, n=n, d=d, p=p, q=q)
#   #   
#   elboPath[iter] <- ElboVecThetaPsi(thetaPsi=c(thetaNew, psiNew), data=data, q=q)
#   diff <- max(abs(c(thetaNew, psiNew) - c(theta, psi)))  
#   theta <- thetaNew; mStep <- Theta2Mstep(theta, n=n, d=d, p=p, q=q)
#   psi <- psiNew; eStep <- Psi2Estep(psi, n=n, d=d, p=p, q=q)
#   if(iter%%round(sqrt(iterMax))==0){
#     cat(iter, ':', diff, fitVE$value, fitM$value, elboPath[iter], vem$elbo, vem$iter, '\n')
#     plot(elboPath[1:iter], type='b', ylim=quantile(elboPath, prob=c(0.1, 1), na.rm=TRUE))
#   }
# }
# elboPath <- elboPath[1:iter]
# plot(elboPath, type='b', ylim=quantile(elboPath, prob=c(0.1, 1)))

# Loop over sims
for(oo in 1:obsNb){
  for(seed in seedList){
    obs <- obsList[[oo]]; set.seed(seed)
    # Data
    simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-coefC', coefC, '-seed', seed)
    simParms <- paste0(simParmsFull, '-obs', 100*obs)
    simName <- paste0(baseSimName, simParms)
    simFile <- paste0(simDir, simName, '.Rdata')
    load(simFile)
    # Fit
    fitName <- paste0(baseFitName, simParms, '-vec')
    fitFile <- paste0(simDir, fitName, '.Rdata')
    if(!file.exists(fitFile)){
      print(simName)
      init <- InitZiPLN(data, q=q)
      vem <- VemZiPLNvec(data, init=init)
      save(init, vem, file=fitFile)
    }
  }
}

