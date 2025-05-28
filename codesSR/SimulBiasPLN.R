# Sim and fit PLN / PLN-PCA

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
seedParm <- 2; set.seed(seedParm)

# seed <- .Random.seed
library(mvtnorm)
library(PLNmodels)
simDir <- '../simBiasPLN/'

# Dims
d <- 2; p <- 3; q <- 5; coefC <- 1
baseSimName <- 'BiasPLN'
seedList <- 1:100; seedNb <- length(seedList)
nList <- c(100, 500, 1000); nNb <- length(nList)

# Sim loop
for(n in nList){ # n <- min(nList)
  parmsName <- paste0('n', n, '-p', p, '-d', d, '-coefC', coefC, '-seed', seedParm)
  parmFile <- paste0(simDir, baseSimName, parmsName, '-parms.Rdata')
  if(!file.exists(parmFile)){
    X0 <- matrix(rnorm(n*d), n, d); X0[, 1] <- 1
    C <- coefC*matrix(rnorm(p*q), p, q)/sqrt(q)
    Sigma <- C%*%t(C)
    beta0 <- 2
    beta <- matrix(rnorm(d*p), d, p)/sqrt(d); beta[, 1] <- beta[, 1] + beta0
    save(X0, beta, C, Sigma, file=parmFile)
  }else{load(parmFile)}
  for(seed in seedList){ # seed <- 1
    set.seed(seed)
    simFile <- paste0(simDir, baseSimName, parmsName, '-sim', seed, '.Rdata')
    if(!file.exists(simFile)){
      X <- X0
      Z <- rmvnorm(n, sigma=Sigma)
      Lambda <- X0%*%beta + Z
      Y <- matrix(rpois(n*p, exp(Lambda)), n, p)
      # plot(exp(Lambda), 1+Y, log='xy'); abline(0, 1)
      save(X, Z, Y, file=simFile)
    }
  }
}

# Fit loop
for(n in nList){ # n <- min(nList)
  cat('n =', n, ': ')
  parmsName <- paste0('n', n, '-p', p, '-d', d, '-coefC', coefC, '-seed', seedParm)
  for(seed in seedList){ # seed <- 1
    cat(seed, '')
    simFile <- paste0(simDir, baseSimName, parmsName, '-sim', seed, '.Rdata')
    plnFile <- paste0(simDir, baseSimName, parmsName, '-sim', seed, '-pln.Rdata')
    pcaFile <- paste0(simDir, baseSimName, parmsName, '-sim', seed, '-pca.Rdata')
    if(!file.exists(pcaFile)){
      load(simFile)
      pln <- PLN(Y ~ - 1 + X, control=PLN_param(trace=0))
      save(pln, file=plnFile)
      pca <- PLNPCA(Y ~ - 1 + X, control=PLN_param(trace=0), ranks=1:p)
      save(pca, file=pcaFile)
    }
  }
  cat('\n')
}

# Res loop
qNorm <- qnorm(p=(1:seedNb)/(seedNb+1))
betaPlnStat <- betaPcaStat <- sigmaPlnMSE <- sigmaPcaMSE <- varPlnStat <- varPcaStat <- list()
for(nn in 1:nNb){ # nn <- 1
  n <- nList[nn]
  parmsName <- paste0('n', n, '-p', p, '-d', d, '-coefC', coefC, '-seed', seedParm)
  parmFile <- paste0(simDir, baseSimName, parmsName, '-parms.Rdata')
  load(parmFile)
  betaPln <- betaPca <- matrix(NA, seedNb, d*p)
  sigmaPln <- sigmaPca <- matrix(NA, seedNb, p^2)
  varPln <- varPca <- matrix(NA, seedNb, p)
  for(seed in seedList){ # seed <- 1
    plnFile <- paste0(simDir, baseSimName, parmsName, '-sim', seed, '-pln.Rdata')
    load(plnFile)
    betaPln[seed, ] <- as.vector(pln$model_par$B)
    sigmaPln[seed, ] <- as.vector(pln$model_par$Sigma)
    varPln[seed, ] <- diag(pln$model_par$Sigma)
    pcaFile <- paste0(simDir, baseSimName, parmsName, '-sim', seed, '-pca.Rdata')
    load(pcaFile)
    betaPca[seed, ] <- as.vector(pca$models[[min(p, q)]]$model_par$B)
    sigmaPca[seed, ] <- as.vector(pca$models[[min(p, q)]]$model_par$Sigma)
    varPca[seed, ] <- diag(pca$models[[min(p, q)]]$model_par$Sigma)
  }
  cat('\n')
  # Beta
  betaTrue <- rep(1, seedNb)%o%as.vector(beta)
  betaPlnSd <- apply(betaPln, 2, sd); betaPcaSd <- apply(betaPca, 2, sd)
  betaPlnStat[[nn]] <- (betaPln - betaTrue)/(rep(1, seedNb)%o%betaPlnSd)
  betaPcaStat[[nn]] <- (betaPca - betaTrue)/(rep(1, seedNb)%o%betaPcaSd)
  # Sigma
  sigmaTrue <- rep(1, seedNb)%o%as.vector(Sigma)
  sigmaPlnMSE[[nn]] <- rowMeans((sigmaPln - sigmaTrue)^2)
  sigmaPcaMSE[[nn]] <- rowMeans((sigmaPca - sigmaTrue)^2)
  # Var
  varTrue <- rep(1, seedNb)%o%diag(Sigma)
  varPlnSd <- apply(varPln, 2, sd); varPcaSd <- apply(varPca, 2, sd)
  varPlnStat[[nn]] <- (varPln - varTrue)/(rep(1, seedNb)%o%varPlnSd)
  varPcaStat[[nn]] <- (varPca - varTrue)/(rep(1, seedNb)%o%varPcaSd)
}

# Plot loop
# par(mfcol=c(p+1, nNb), mex=.5)
par(mfrow=c(2*nNb, p), mex=.5, pch=20)
for(nn in 1:nNb){ # nn <- 1
  # Beta
  for(j in 1:p){ # j <- 1
    plot(qNorm, sort(betaPlnStat[[nn]][, d*(j-1)+1]), xlab='', ylab='', main=paste0(nList[nn], ' beta1', j))
    points(qNorm, sort(betaPcaStat[[nn]][, d*(j-1)+1]), col=2)
    abline(a=0, b=1, h=0, v=0)
  }
  # # Sigma
  # boxplot(cbind(sigmaPlnMSE[[nn]], sigmaPcaMSE[[nn]]), log='y', xlab='', ylab='', main=paste0(nList[nn], ' sigmaMS'))
  # Var
  for(j in 1:p){ # j <- 1
    plot(qNorm, sort(varPlnStat[[nn]][, j]), xlab='', ylab='', main=paste0(nList[nn], ' sigma', j, j))
    points(qNorm, sort(varPcaStat[[nn]][, j]), col=2)
    abline(a=0, b=1, h=0, v=0)
  }
}

par(mfcol=c(p, nNb*(d-1)), mex=.5, pch=20)
for(nn in 1:nNb){ # nn <- 1
  # Beta
  for(h in 2:d){
    for(j in 1:p){ # j <- 1
      plot(qNorm, sort(betaPlnStat[[nn]][, d*(j-1)+h]), xlab='', ylab='', main=paste0(nList[nn], ' beta', h, j))
      points(qNorm, sort(betaPcaStat[[nn]][, d*(j-1)+h]), col=2)
      abline(a=0, b=1, h=0, v=0)
    }
  }
}
