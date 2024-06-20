# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
library(Rcpp)
library(PLNmodels)

# seed <- .Random.seed
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIPLN.R')
source('codesSR/Functions/FunctionsZIP.R')
source("FunctionsBB.R")
simDir <- 'SimulationsBB/datasim/'
# Parms: many small sims
n <- 100; d <- 5; p <- 10; q <- 2
baseSimName <- 'ZiPLNsim'
seedList <- 1:10; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)

# # Parms: one big sim
# n <- 500; d <- 20; p <- 30; q <- 5
# seedList <- 1; seedNb <- length(seedList)
# obsList <- c(0.6); obsNb <- length(obsList)

#X0 <- NULL
 # Same X for all
seedX <- 1; set.seed(seedX)
#X0 <- matrix(rnorm(n*p*d), n*p, d); X0[, 1] <- 1; 
baseSimName <- paste0(baseSimName, '-sameX', seedX)

# Simul
for(seed in 1:10){
  set.seed(seed)
  simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
  simNameFull <- paste0(baseSimName, simParmsFull)
  simFileFull <- paste0(simDir, simNameFull, '-noMiss.Rdata')
  if(!file.exists(simFileFull)){
    sim <- SimZiPLN(n=n, p=p, d=d, q=q, X=X0)
    obsTresh <- matrix(runif(n*p), n, p)
    data <- list(X=sim$X, Y=sim$Y, obsTresh=obsTresh, ij=sim$ij, logFactY=lgamma(sim$Y+1))
    true <- list(mStep=list(gamma=sim$gamma, beta=sim$beta, C=sim$C),
                 eStep=list(M=sim$W, S=matrix(1e-4, n, q), xi=matrix(plogis(sim$X%*%sim$gamma), n, p)),
                 latent=list(U=sim$U, W=sim$W, Z=sim$Z, Yall=sim$Yall))
    save(data, true, file=simFileFull)
    for(oo in 1:obsNb){
      obs <- obsList[oo]
      simParms <- paste0(simParmsFull, '-obs', 100*obs)
      simName <- paste0(baseSimName, simParms)
      simFile <- paste0(simDir, simName, '.Rdata')
      Omega <- 1*(obsTresh <= obs)
      data <- list(X=sim$X, Y=sim$Y, Omega=Omega, ij=sim$ij, logFactY=lgamma(sim$Y+1))
      save(data, file=simFile)
    }
  }
}




# Fit de Stephane
start_time <- Sys.time()
for(seed in 1:1){
    simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
    simNameFull <- paste0(baseSimName, simParmsFull)
    simFileFull <- paste0(simDir, simNameFull, '-noMiss.Rdata')
    load(file=simFileFull)
    
    
    ################ ajustement SR
    init <- InitZiPLN(data) #! Fichiers différnets pour miss et pas miss
    vem <- VemZiPLN(data=data, init=init, iterMax=5e3)
    # save(init, vem, file=fitFile)
    ################ ajustement BB
    
    for(oo in 1:obsNb){
        obs <- obsList[oo]
        simParms <- paste0(simParmsFull, '-obs', 100*obs)
        simName <- paste0(baseSimName, simParms)
        simFile <- paste0(simDir, simName, '.Rdata')
      ################ ajustement SR
       
        
    }
    }

end_time <- Sys.time()
execution_time_real <- end_time - start_time

#############################################################
# Fit de Barbara

params <- list(B = as.matrix(init$mStep$beta), 
               D = as.matrix(init$mStep$gamma), 
               C = as.matrix(init$mStep$C), 
               M = as.matrix(init$eStep$M), 
               S = as.matrix(init$eStep$S))


source("FunctionsBB.R")

config <- PLNPCA_param()$config_optim

# start_timeB <- Sys.time()
# outB <- Miss.ZIPPCA(data$Y, data$X, q, params = params)
# end_timeB <- Sys.time()
# execution_time_realB <- end_timeB - start_timeB
# 
# outB$elbo
# vem$elbo
# 
# plot(log(1 + data$Y[data$Y != 0]), log(1 + vem$pred$A[data$Y != 0])) ; abline(0,1)
# plot(log(1 + data$Y[data$Y != 0]), log(1 + outB$pred$A[data$Y != 0])) ; abline(0,1)
# plot(log(1 + vem$pred$A[data$Y != 0]), log(1 + outB$pred$A[data$Y != 0])) ; abline(0,1)
# 
# 
# boxplot(vem$pred$nu[data$Y == 0], vem$pred$nu[data$Y != 0])
# 
# boxplot(outB$pred$nu[data$Y == 0], outB$pred$nu[data$Y != 0])

data$Omega <- ifelse(is.na(data$Y), 0, 1)
data$R <- ifelse(is.na(data$Y), 0, 1)

source("codesSR/Functions/FunctionsZIPLNmiss.R")
source("codesSR/Functions/FunctionsZIPLNmissVec.R")


ElboSfun <- ELBO(data = data, mStep = init$mStep, eStep = init$eStep)
SgradS <- matrix(ElboGradVecS(Svec = as.vector(t(init$eStep$S)), data = data, mStep = init$mStep, eStep = init$eStep),n, q, byrow = TRUE)
SgradM <- matrix(ElboGradVecM(Mvec = as.vector(t(init$eStep$M)), data = data, mStep = init$mStep, eStep = init$eStep),n, q, byrow = TRUE)
SgradBeta <- ElboGradBeta(beta = init$mStep$beta, data = data, mStep = init$mStep, eStep = init$eStep)
SgradGamma <- ElboGradGamma(gamma = init$mStep$gamma, data = data, mStep = init$mStep, eStep = init$eStep)
SgradC <- as.matrix(ElboGradC(vecC = as.vector(init$mStep$C), data = data, mStep = init$mStep, eStep = init$eStep),p,q)

ElboBfun.init <- ElboB(data = data, params = params)

BgradS <- ElboBfun.init$gradS
BgradM <-ElboBfun.init$gradM
BgradBeta <- ElboBfun.init$gradB
BgradGamma <- ElboBfun.init$gradD
BgradC <- ElboBfun.init$gradC

BgradS - SgradS
BgradM - SgradM












