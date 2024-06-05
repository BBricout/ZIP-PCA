# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

# seed <- .Random.seed
source('Functions/FunctionsUtils.R')
source('Functions/FunctionsZIP.R')
source('Functions/FunctionsZIPLNmiss.R')
source('Functions/FunctionsZIPLNmissVec.R')
simDir <- '../simulSR/'

# Parms: 
n <- 100; d <- 5; p <- 10; q <- 2; obs <- 0.70; seed <- 1
simParms <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed, '-obs', 100*obs)
simName <- 'ZiPLNsim'; fitName <- 'ZiPLNfit'; 

# 1 example
load(paste0(simDir, simName, simParms, '.Rdata'))
load(paste0(simDir, fitName, simParms, '.Rdata'))
vem$eStep$xi <- ComputeXi(data=data, mStep=vem$mStep, eStep=vem$eStep)
vem$elbo <- ELBO(data=data, mStep=vem$mStep, eStep=vem$eStep)

# # Init
# initThetaPsi <- as.vector(c(init$mStep$gamma, init$mStep$beta, as.vector(init$mStep$C), 
#                             as.vector(t(init$eStep$M)), as.vector(t(init$eStep$S))))
# ElboVecThetaPsi(thetaPsi=initThetaPsi, data=data, q=q)
# ElboGradVecThetaPsi(thetaPsi=initThetaPsi, data=data, q=q)

# # Vem
# vemThetaPsi <- as.vector(c(vem$mStep$gamma, vem$mStep$beta, as.vector(vem$mStep$C), 
#                   as.vector(t(vem$eStep$M)), as.vector(t(vem$eStep$S))))
# ElboVecThetaPsi(thetaPsi=vemThetaPsi, data=data, q=q)
# ElboGradVecThetaPsi(thetaPsi=vemThetaPsi, data=data, q=q)

# VEM
mStep <- init$mStep; theta <- Mstep2Theta(mStep, n=n, d=d, p=p, q=q)
eStep <- init$eStep; psi <- Estep2Psi(eStep, n=n, d=d, p=p, q=q)
tolS <- 1e-4; iterMax <- 1e2; tol <- 1e-4; tolXi <- 1e-4
elboPath <- rep(NA, iterMax)
iter <- 1; diff <- 2*tol
elboPath[iter] <- ElboVecThetaPsi(thetaPsi=c(theta, psi), data=data, q=q)
while((iter & iterMax) & (diff > tol)){
  iter <- iter+1
  # # VE + M
  # fit <- optim(par=c(theta, psi), fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=data,
  #              q=q, method='L-BFGS-B', control=list(fnscale=-1),
  #              lower=c(rep(-Inf, ((2*d) + (p*q) + (n*q))), rep(tolS, (n*q))))
  # thetaNew <- fit$par[1:((2*d)+(p*q))]; psiNew <- fit$par[-(1:((2*d)+(p*q)))]
  # VE
  fitVE <- optim(par=psi, fn=ElboVecPsi, gr=ElboGradVecPsi, data=data, mStep=mStep,
               method='L-BFGS-B', control=list(fnscale=-1),
               lower=c(rep(-Inf, (n*q)), rep(tolS, (n*q))))
  psiNew <- fitVE$par; eStepNew <- Psi2Estep(psiNew, n=n, d=d, p=p, q=q)
  eStepNew$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStepNew, tolXi=tolXi)
  # M
  fitM <- optim(par=theta, fn=ElboVecTheta, gr=ElboGradVecTheta, data=data, eStep=eStepNew, 
               method='BFGS', control=list(fnscale=-1))
  thetaNew <- fitM$par; mStepNew <- Theta2Mstep(thetaNew, n=n, d=d, p=p, q=q)
  #   
  elboPath[iter] <- ElboVecThetaPsi(thetaPsi=c(thetaNew, psiNew), data=data, q=q)
  diff <- max(abs(c(thetaNew, psiNew) - c(theta, psi)))  
  theta <- thetaNew; mStep <- Theta2Mstep(theta, n=n, d=d, p=p, q=q)
  psi <- psiNew; eStep <- Psi2Estep(psi, n=n, d=d, p=p, q=q)
  cat(iter, ':', diff, fitVE$value, fitM$value, elboPath[iter], vem$elbo, vem$iter, '\n')
}
elboPath <- elboPath[1:iter]
plot(elboPath, type='b', ylim=quantile(elboPath, prob=c(0.1, 1)))
