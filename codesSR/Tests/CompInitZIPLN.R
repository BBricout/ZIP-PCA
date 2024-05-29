# Compares init for ZIPLN

# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
# seed <- .Random.seed
source('../FunctionsZIP.R'); 
source('../FunctionsZIPLN.R')
simDir <- '../../simulSR/'

# Parms
n <- 100; d <- 5; p <- 10; q <- 2; seed <- 8
simParms <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
simName <- paste0('ZiPLNsim', simParms)
simFile <- paste0(simDir, simName, '.Rdata')
load(simFile)
data <- list(X=sim$X, Y=sim$Y, ij=sim$ij, logFactY=lgamma(sim$Y+1))

# Inits
init0 <- InitZiPLNold(data)
init1 <- InitZiPLN(data)

# Compares
c(ELBO(data=data, mStep=init0$mStep, eStep=init0$eStep), 
  ELBO(data=data, mStep=init1$mStep, eStep=init1$eStep))

rbind(sim$gamma, init0$mStep$gamma, init1$mStep$gamma)
plot(as.data.frame(cbind(data$X%*%sim$gamma, data$X%*%init0$mStep$gamma, data$X%*%init1$mStep$gamma)))
rbind(sim$beta, init0$mStep$beta, init1$mStep$beta)
plot(as.data.frame(cbind(data$X%*%sim$beta, data$X%*%init0$mStep$beta, data$X%*%init1$mStep$beta)))
