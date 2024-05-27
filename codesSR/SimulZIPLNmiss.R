# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

# seed <- .Random.seed
source('Functions/FunctionsUtils.R')
source('Functions/FunctionsZIPLNmiss.R')
simDir <- '../simulSR/'

# Parms: many small sims
n <- 100; d <- 5; p <- 10; q <- 2
seedList <- 1:10; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)

# # Parms: one big sim
# n <- 500; d <- 20; p <- 30; q <- 5
# seedList <- 1; seedNb <- length(seedList)
# obsList <- c(0.6); obsNb <- length(obsList)

# Simul
for(seed in 1:10){
  set.seed(seed)
  simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
  simNameFull <- paste0('ZiPLNsim', simParmsFull)
  simFileFull <- paste0(simDir, simNameFull, '-noMiss.Rdata')
  if(!file.exists(simFileFull)){
    sim <- SimZiPLN(n=n, p=p, d=d, q=q)
    obsTresh <- matrix(runif(n*p), n, p)
    data <- list(X=sim$X, Y=sim$Y, obsTresh=obsTresh, ij=sim$ij, logFactY=lgamma(sim$Y+1))
    true <- list(mstep=list(gamma=sim$gamma, beta=sim$beta, C=sim$C), 
                 eStep=list(M=sim$W, S=matrix(1e-4, n, q), xi=matrix(plogis(sim$X%*%sim$gamma), n, p)), 
                 latent=list(U=sim$U, W=sim$W, Z=sim$Z, Yall=sim$Yall))
    save(data, true, file=simFileFull)
    for(oo in 1:obsNb){
      obs <- obsList[oo]
      simParms <- paste0(simParmsFull, '-obs', 100*obs)
      simName <- paste0('ZiPLNsim', simParms)
      simFile <- paste0(simDir, simName, '.Rdata')
      Omega <- 1*(obsTresh <= obs)
      data <- list(X=sim$X, Y=sim$Y, Omega=Omega, ij=sim$ij, logFactY=lgamma(sim$Y+1))
      save(data, file=simFile)
    }
  }
}
