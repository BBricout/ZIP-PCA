# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

# seed <- .Random.seed
source('Functions/FunctionsUtils.R')
source('Functions/FunctionsZIPLNmiss.R')
simDir <- '../simulSR/'

# Parms: many small sims
n <- 100; d <- 2; p <- 5; q <- 2; coefC <- 1
# n <- 100; d <- 5; p <- 10; q <- 2; obs <- 0.6
baseSimName <- 'ZiPLNsim'; baseFitName <- 'ZiPLNfit'; 
seedList <- 1:100; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)
obsList <- c(1); obsNb <- length(obsList)

# Functions
MakeCOrtho <- function(C){
  q <- ncol(C); eigC <- eigen(C%*%t(C))
  return(eigC$vectors[, 1:q]%*%diag(sqrt(eigC$values[1:q])))
}

# # Parms: one big sim
# n <- 500; d <- 20; p <- 30; q <- 5
# seedList <- 1; seedNb <- length(seedList)
# obsList <- c(0.6); obsNb <- length(obsList)

# # Parms
# X0 <- matrix(rnorm(n*p*d), n*p, d); X0[, 1] <- 1
# gamma <- rnorm(d)/sqrt(d)
# C <- coefC*matrix(rnorm(p*q), p, q)/sqrt(q)
# beta0 <- 2
# beta <- rnorm(d)/sqrt(d); beta[1] <- beta[1] + beta0

# Simul
simParms <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-coefC', coefC)
thetaTrue <- thetaHat <- thetaSim <- matrix(NA, 100, 2*d)
colnames(thetaHat) <- colnames(thetaTrue) <- colnames(thetaSim) <- c(paste0('gamma', 1:d), paste0('beta', 1:d))
for(seed in seedList){ # seed <- 7
  parmsName <- paste0(simParms, '-seed', seed)
  simFileFull <- paste0(simDir, baseSimName, parmsName, '-noMiss.Rdata')
  if(file.exists(simFileFull)){
    load(simFileFull)
    thetaTrue[seed, ] <- c(true$mStep$gamma, true$mStep$beta)
    for(oo in 1:1){ # oo <- 1
      obs <- obsList[oo]
      simName <- paste0(parmsName, '-obs', 100*obs)
      fitFile <- paste0(simDir, baseFitName, simName, '.Rdata')
      load(fitFile)
      thetaHat[seed, ] <- c(vem$mStep$gamma, vem$mStep$beta)
      vemOrth <- vem; vemOrth$mStep$C <- MakeCOrtho(vem$mStep$C)
      ELBO(data, true$mStep, true$eStep)
      ELBO(data, vem$mStep, vem$eStep)
      ELBO(data, vemOrth$mStep, vemOrth$eStep)
      nuMuA <- NuMuA(data, vem$mStep, vem$eStep)
      nuMuAorth <- NuMuA(data, vemOrth$mStep, vemOrth$eStep)
      plot(nuMuA$A, nuMuAorth$A)
    }
  }
}

boxplot(thetaHat-thetaTrue); abline(h=0, col=2)
thetaSd <- apply(thetaHat, 2, sd)
par(mfrow=c(2, d))
for(h in 1:(2*d)){
  qqnorm((thetaHat[, h] - thetaTrue[, h])/thetaSd[h])
  abline(a=0, b=1, h=0, v=0)
}
