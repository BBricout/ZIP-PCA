# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

# seed <- .Random.seed
source('Functions/FunctionsUtils.R')
source('Functions/FunctionsZIPLNmiss.R')
simDir <- '../simulSR/'

# Parms: many small sims
n <- 500; d <- 2; p <- 5; q <- 2; coefC <- 1
# n <- 100; d <- 5; p <- 10; q <- 2; obs <- 0.6
baseSimName <- 'ZiPLNsim'
seedList <- 1:100; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)

# # Parms: one big sim
# n <- 500; d <- 20; p <- 30; q <- 5
# seedList <- 1; seedNb <- length(seedList)
# obsList <- c(0.6); obsNb <- length(obsList)

# Parms
X0 <- matrix(rnorm(n*p*d), n*p, d); X0[, 1] <- 1
gamma <- rnorm(d)/sqrt(d)
C <- coefC*matrix(rnorm(p*q), p, q)/sqrt(q)
beta0 <- 2
beta <- rnorm(d)/sqrt(d); beta[1] <- beta[1] + beta0

# Simul
simParms <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-coefC', coefC)
for(seed in 1:100){ # seed <- 1
  set.seed(seed)
  parmsName <- paste0(baseSimName, simParms, '-seed', seed)
  simFileFull <- paste0(simDir, parmsName, '-noMiss.Rdata')
  if(!file.exists(simFileFull)){
    sim <- SimZiPLN(n=n, p=p, d=d, q=q, X=X0, coefC=coefC, gamma=gamma, beta=beta, C=C)
    obsTresh <- matrix(runif(n*p), n, p)
    data <- list(X=sim$X, Y=sim$Y, obsTresh=obsTresh, ij=sim$ij, logFactY=lgamma(sim$Y+1))
    true <- list(mStep=list(gamma=sim$gamma, beta=sim$beta, C=sim$C), 
                 eStep=list(M=sim$W, S=matrix(1e-4, n, q), xi=matrix(plogis(sim$X%*%sim$gamma), n, p)), 
                 latent=list(U=sim$U, W=sim$W, Z=sim$Z, Yall=sim$Yall))
    save(data, true, file=simFileFull)
    for(oo in 1:obsNb){ # oo <- 1
      obs <- obsList[oo]
      simName <- paste0(parmsName, '-obs', 100*obs)
      simFile <- paste0(simDir, simName, '.Rdata')
      Omega <- 1*(obsTresh <= obs)
      data <- list(X=sim$X, Y=sim$Y, Omega=Omega, ij=sim$ij, logFactY=lgamma(sim$Y+1))
      save(data, file=simFile)
    }
  }
}

