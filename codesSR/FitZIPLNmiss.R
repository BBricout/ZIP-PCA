# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

# seed <- .Random.seed
source('Functions/FunctionsUtils.R')
source('Functions/FunctionsZIP.R')
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

# Loop over sims
for(seed in seedList){
  for(oo in 1:obsNb){
    obs <- obsList[[oo]]; set.seed(seed)
    # Data
    simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
    simParms <- paste0(simParmsFull, '-obs', 100*obs)
    simName <- paste0('ZiPLNsim', simParms)
    simFile <- paste0(simDir, simName, '.Rdata')
    load(simFile)
    # Fit
    fitName <- paste0('ZiPLNfit', simParms)
    fitFile <- paste0(simDir, fitName, '.Rdata')
    if(!file.exists(fitFile)){
      print(simName)
      init <- InitZiPLN(data)
      vem <- VemZiPLN(data, init=init)
      save(init, vem, file=fitFile)
    }
  }
}
