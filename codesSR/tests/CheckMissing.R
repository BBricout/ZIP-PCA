# Compares init for ZIPLN

# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
# seed <- .Random.seed
source('../FunctionsZIP.R'); 
source('../FunctionsZIPLNmiss.R')
simDir <- '../../simulSR/'

# Parms
n <- 100; d <- 5; p <- 10; q <- 2; seed <- 2; obs <- 0.99
simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
simParms <- paste0(simParmsFull, '-obs', 100*obs)
simName <- paste0('ZiPLNsim', simParms)
simFile <- paste0(simDir, simName, '.Rdata')
load(simFile)

# Compare init
data0 <- list(X=sim$X, Y=sim$Y, Omega=sim$Omega, ij=sim$ij, logFactY=lgamma(sim$Y+1))
init0 <- InitZiPLN(data=data0)
elbo0 <- ELBO(data=data0, mStep=init0$mStep, eStep=init0$eStep)

sim$Y[which(sim$Omega==0)] <- 12
data1 <- list(X=sim$X, Y=sim$Y, Omega=sim$Omega, ij=sim$ij, logFactY=lgamma(sim$Y+1))
init1 <- InitZiPLN(data=data1)
elbo1 <- ELBO(data=data1, mStep=init1$mStep, eStep=init1$eStep)

rbind(init0$mStep$gamma, init1$mStep$gamma)
max(abs(init0$mStep$gamma - init1$mStep$gamma))
rbind(init0$mStep$beta, init1$mStep$beta)
max(abs(init0$mStep$beta - init1$mStep$beta))
c(elbo0, elbo1)

# Compare Mstep
mStep0 <- Mstep(data=data0, mStep=init0$mStep, eStep=init0$eStep)
mStep1 <- Mstep(data=data1, mStep=init1$mStep, eStep=init1$eStep)
rbind(mStep0$gamma, mStep0$gamma)
max(abs(mStep0$gamma - mStep1$gamma))
rbind(mStep1$beta, mStep1$beta)
max(abs(mStep0$beta - mStep1$beta))

# Compare Estep
par(mfrow=c(2, 2))
eStep0 <- VEstep(data=data0, mStep=init0$mStep, eStep=init0$eStep)
eStep1 <- VEstep(data=data1, mStep=init1$mStep, eStep=init1$eStep)
plot(mStep0$C, mStep1$C, main='C'); abline(0, 1)
max(abs(mStep0$C - mStep1$C))
plot(eStep0$xi, eStep1$xi, main='xi'); abline(0, 1)
max(abs(eStep0$xi - eStep1$xi))
plot(eStep0$M, eStep1$M, main='M'); abline(0, 1)
max(abs(eStep0$M - eStep1$M))
plot(eStep0$S, eStep1$S, main='S'); abline(0, 1)
max(abs(eStep0$S - eStep1$S))
