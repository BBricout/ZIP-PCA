# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

# seed <- .Random.seed
source('Functions/FunctionsUtils.R')
source('Functions/FunctionsZIP.R')
source('Functions/FunctionsZIPLNmiss.R')
source('Functions/FunctionsZIPLNmissVec.R')
simDir <- '../simulSR/'

# Parms: many small sims
n <- 100; d <- 5; p <- 10; q <- 2
baseSimName <- 'ZiPLNsim'; baseFitName <- 'ZiPLNfit'; 
# baseSimName <- 'ZiPLNsim-sameX1'; baseFitName <- 'ZiPLNfit-sameX1';
seedList <- 1:10; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)

# # Parms: one big sim
# n <- 500; d <- 20; p <- 30; q <- 5
# seedList <- 1; seedNb <- length(seedList)
# obsList <- c(0.6); obsNb <- length(obsList)

# 1 example
simFile <- "ZiPLNsim-n100-d5-p10-q2-seed1-obs95.Rdata"
load(paste0(simDir, simFile))
fitFile <- "ZiPLNfit-n100-d5-p10-q2-seed1-obs95.Rdata"
load(paste0(simDir, fitFile))
vem$eStep$xi <- ComputeXi(data=data, mStep=vem$mStep, eStep=vem$eStep)
vem$elbo <- ELBO(data=data, mStep=vem$mStep, eStep=vem$eStep)

initThetaPsi <- as.vector(c(init$mStep$gamma, init$mStep$beta, as.vector(init$mStep$C), 
                            as.vector(t(init$eStep$M)), as.vector(t(init$eStep$S))))
initTheta <- c(init$mStep$gamma, init$mStep$beta, as.vector(init$mStep$C))
initPsi <- c(as.vector(t(init$eStep$M)), as.vector(t(init$eStep$S)))
eStep <- init$eStep; mStep <- init$mStep
ElboVecThetaPsi(thetaPsi=initThetaPsi, data=data, q=q)
ElboGradVecThetaPsi(thetaPsi=initThetaPsi, data=data, q=q)

# vemThetaPsi <- as.vector(c(vem$mStep$gamma, vem$mStep$beta, as.vector(vem$mStep$C), 
#                   as.vector(t(vem$eStep$M)), as.vector(t(vem$eStep$S))))
# ElboVecThetaPsi(thetaPsi=vemThetaPsi, data=data, q=q)
# ElboGradVecThetaPsi(thetaPsi=vemThetaPsi, data=data, q=q)

tolS <- 1e-4
ElboVecThetaPsi(thetaPsi=initThetaPsi, data=data, q=q)
# fit <- optim(par=initThetaPsi, fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=data,
#              q=q, method='L-BFGS-B', control=list(fnscale=-1),
#              lower=c(rep(-Inf, ((2*d) + (p*q) + (n*q))), rep(tolS, (n*q))))
# fit$value
fit <- optim(par=initPsi, fn=ElboVecPsi, gr=ElboGradVecPsi, data=data, mStep=mStep,
             method='L-BFGS-B', control=list(fnscale=-1), 
             lower=c(rep(-Inf, (n*q)), rep(tolS, (n*q))))
fit$value; 
eStep <- Psi2Estep(fit$par, n=n, d=d, p=p, q=q)
eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=0)
fit <- optim(par=initTheta, fn=ElboVecTheta, gr=ElboGradVecTheta, data=data, eStep=eStep, 
             method='BFGS', control=list(fnscale=-1))
fit$value
mStep <- Theta2Mstep(fit$par, n=n, d=d, p=p, q=q)

# Loop over sims
for(seed in seedList){
  for(oo in 1:obsNb){
    obs <- obsList[[oo]]; set.seed(seed)
    # Data
    simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
    simParms <- paste0(simParmsFull, '-obs', 100*obs)
    simName <- paste0(baseSimName, simParms)
    simFile <- paste0(simDir, simName, '.Rdata')
    load(simFile)
    # Fit
    fitName <- paste0(baseFitName, simParms)
    fitFile <- paste0(simDir, fitName, '.Rdata')
    if(!file.exists(fitFile)){
      print(simName)
      init <- InitZiPLN(data)
      vem <- VemZiPLN(data, init=init)
      save(init, vem, file=fitFile)
    }
  }
}

