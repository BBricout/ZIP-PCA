# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

# seed <- .Random.seed
source('FunctionsZIP.R')
source('FunctionsZIPLNmiss.R')
simDir <- '../simulSR/'

# Parms
n <- 100; d <- 5; p <- 10; q <- 2
seedList <- 1:10; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)
for(seed in seedList){
  for(oo in 1:obsNb){
    obs <- obsList[[oo]]; set.seed(seed)
    # Data
    simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
    simNameFull <- paste0('ZiPLNsim', simParmsFull)
    simFileFull <- paste0(simDir, simNameFull, '.Rdata')
    load(simFileFull)
    simParms <- paste0(simParmsFull, '-obs', 100*obs)
    simName <- paste0('ZiPLNsim', simParms)
    simFile <- paste0(simDir, simName, '.Rdata')
    if(!file.exists(simFile)){
      sim <- SimZiPLNmiss(sim=sim, obs=obs)
      save(sim, file=simFile)
    }else{load(simFile)}
    data <- list(X=sim$X, Y=sim$Y, Omega=sim$Omega, ij=sim$ij, logFactY=lgamma(sim$Y+1))
    # Fit
    fitName <- paste0('ZiPLNfit', simParms)
    fitFile <- paste0(simDir, fitName, '.Rdata')
    if(!file.exists(fitFile)){
      print(simName)
      oracle <- OracleZiPLN(sim)
      init <- InitZiPLN(data)
      vem <- VemZiPLN(data, init=init)
      save(oracle, init, vem, file=fitFile)
    }
  }
}
