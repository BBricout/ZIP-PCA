################################################################################
# Criteria
Criteria <- function(data, fit){
  n <- nrow(data$Y); d <- ncol(data$X); p <- ncol(data$Y); q <- ncol(fit$eStep$M)
  penBic <- (2*d + choose(p, 2) - choose(q-1, 2))*log(n)/2
  ent <- 0.5*(n*q + log(2*pi*exp(1)) + sum(log(fit$eStep$S)))
  ent <- ent - sum(fit$eStep$xi*log(fit$eStep$xi + (fit$eStep$xi==0)) + 
                     (1 - fit$eStep$xi)*log(1 - fit$eStep$xi + (fit$eStep$xi==1)))
  if(is.null(fit$diff)){fit$diff <- NA}
  return(list(iter=fit$iter, diff=fit$diff, elbo=fit$elbo, penBic=penBic, bic=fit$elbo-penBic, 
              ent=ent, icl=fit$elbo-penBic-ent))
}


################################################################################
# Jackknife
JackknifeZiPLN <- function(data, fit, iterMax=1e3){
  n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X); q <- ncol(fit$eStep$M)
  fit_iList <- list()
  ij_1 <- cbind(rep(1:(n-1), p), rep(1:p, each=n-1))
  gammaMat <- betaMat <- matrix(NA, n, d)
  CMat <- matrix(NA, n, p*q)
  xiMat <- matrix(NA, n, p); MMat <- SMat <- matrix(NA, n, q)
  par(mfrow=c(5, 5))
  for(i in 1:n){ # i <- 1
    cat('i =', i, ': ')
    data_i <- list(X=data$X[-which(data$ij[, 1]==i), ], Y=data$Y[-i, ], Omega=data$Omega[-i, ],
                   ij=ij_1, logFactY=data$logFactY[-i, ])
    init_i <- list(mStep=list(gamma=fit$mStep$gamma, beta=fit$mStep$beta, C=fit$mStep$C),
                   eStep=list(xi=fit$eStep$xi[-i, ], M=fit$eStep$M[-i, ], S=fit$eStep$S[-i, ]))
    fit_iList[[i]] <- VemZiPLN(data=data_i, init=init_i, plot=FALSE, iterMax=iterMax)
    plot(fit_iList[[i]]$elboPath, main=i, ylim=quantile(fit_iList[[i]]$elboPath, probs=c(0.1, 1)), xlab='', ylab='')
    gammaMat[i, ] <- fit_iList[[i]]$mStep$gamma; betaMat[i, ] <- fit_iList[[i]]$mStep$beta
    CMat[i, ] <- as.vector(fit_iList[[i]]$mStep$C)
    eStep_i <- VEstep(data=data, mStep=fit_iList[[i]]$mStep, eStep=fit$eStep)
    xiMat[i, ] <- eStep_i$xi[i, ]; MMat[i, ] <- eStep_i$M[i, ]; SMat[i, ] <- eStep_i$S[i, ]
  }
  jkCoef <- (n-1)^2/n
  return(list(fit_iList=fit_iList, 
              mStepMat=list(gamma=gammaMat, beta=betaMat, C=CMat),
              cov=list(gamma=jkCoef*cov(gammaMat), beta=jkCoef*cov(betaMat), 
                       C=jkCoef*cov(CMat), theta=jkCoef*cov(cbind(gammaMat, betaMat, CMat))), 
              eStep=list(xi=xiMat, M=MMat, S=SMat)
              ))
}

################################################################################
# Prediction
PredZiPLN <- function(data, fit){
  n <- nrow(data$Y); p <- ncol(data$Y)
  sigma <- fit$mStep$C %*% t(fit$mStep$C)
  nuMuA <-NuMuA(data=data, mStep=fit$mStep, eStep=fit$eStep)
  condNu <- margNu <- nuMuA$nu
  condPi <- margPi <- plogis(condNu)
  margLambda <- exp(nuMuA$mu + rep(1, n)%o%diag(sigma)/2)
  condLambda <- nuMuA$A
  margEsp <- margPi*margLambda
  margVar <- margEsp + margPi*(1-margPi)*margLambda^2
  condEsp <- condPi*condLambda
  condVar <- condEsp + condPi*(1-condPi)*condLambda^2
  margZipCIl <- margZipCIu <- condZipCIl <- condZipCIu <- 
    margPoisCIl <- margPoisCIu <- condPoisCIl <- condPoisCIu <- matrix(NA, n, p)
  for(i in 1:n){for(j in 1:p){
    margZipCIl[i, j] <- qzip(p=.025, psi=1-margPi[i, j], lambda=margLambda[i, j])
    margZipCIu[i, j] <- qzip(p=.975, psi=1-margPi[i, j], lambda=margLambda[i, j])
    margPoisCIl[i, j] <- qpois(p=.025, lambda=margLambda[i, j])
    margPoisCIu[i, j] <- qpois(p=.975, lambda=margLambda[i, j])
    condZipCIl[i, j] <- qzip(p=.025, psi=1-condPi[i, j], lambda=condLambda[i, j])
    condZipCIu[i, j] <- qzip(p=.975, psi=1-condPi[i, j], lambda=condLambda[i, j])
    condPoisCIl[i, j] <- qpois(p=.025, lambda=condLambda[i, j])
    condPoisCIu[i, j] <- qpois(p=.975, lambda=condLambda[i, j])
  }}
  return(list(nuMuA=nuMuA,
              marg=list(pi=margPi, lambda=margLambda, esp=margEsp, var=margVar, 
                        lZipCI=margZipCIl, uZipCI=margZipCIu, lPoisCI=margPoisCIl, uPoisCI=margPoisCIu), 
              cond=list(pi=condPi, lambda=condLambda, esp=condEsp, var=condVar, 
                        lZipCI=condZipCIl, uZipCI=condZipCIu), lPoisCI=condPoisCIl, uPoisCI=condPoisCIu))
}

