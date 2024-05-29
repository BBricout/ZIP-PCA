# Functions ZI-PLN

################################################################################
# Utils
NuMuA <- function(data, mStep, eStep){
  nu <- matrix(data$X%*%mStep$gamma, nrow(data$Y), ncol(data$Y))
  mu <- matrix(data$X%*%mStep$beta, nrow(data$Y), ncol(data$Y))
  A <- exp(mu + eStep$M%*%t(mStep$C) + 
             0.5*t(sapply(1:nrow(data$Y), function(i){diag(mStep$C%*%diag(eStep$S[i, ])%*%t(mStep$C))})))
  return(list(nu=nu, mu=mu, A=A))
}

################################################################################
# Init
InitZiPLNold <- function(data, q){
  pres <- which(data$Omega==1)
  reg <- lm(as.vector(log(1+data$Y)[pres]) ~ -1 + data$X[pres, ])
  res <- matrix(0, nrow(data$Y), ncol(data$Y))
  res[pres] <- reg$residuals
  pca <- prcomp(res, rank=q)
  mStep <- list(gamma=as.vector(glm(as.vector(1*(data$Y > 0)) ~ -1 + data$X, family='binomial')$coef), 
                beta=reg$coef, 
                C=pca$rotation %*% diag(pca$sdev[1:q]))
  eStep <- list(xi=matrix(sum(data$Omega*(data$Y > 0))/sum(data$Omega), n, ncol(data$Y)), M=matrix(0, n, q), S=matrix(1e-4, n, q))
  return(list(mStep=mStep, eStep=eStep, reg=reg, pca=pca))
}
InitZiPLN <- function(data, q){
  pres <- which(data$Omega==1)
  reg <- lm(as.vector(log(1+data$Y))[pres] ~ -1 + data$X[pres, ])
  res <- matrix(0, nrow(data$Y), ncol(data$Y))
  res[pres] <- reg$residuals
  pca <- prcomp(res, rank=q)
  zip <- EmZIP(X=data$X[pres, ], Y=as.vector(data$Y)[pres])
  mStep <- list(gamma=zip$gamma, beta=zip$beta, 
                C=pca$rotation %*% diag(pca$sdev[1:q]))
  # eStep <- list(xi=matrix(sum(data$Omega*(data$Y > 0))/sum(data$Omega), n, p), 
                # M=matrix(0, n, q), S=matrix(1e-4, n, q))
  xi <- matrix(plogis(data$X%*%mStep$gamma), nrow(data$Y), ncol(data$Y))
  xi[which(data$Y > 0)] <- 1
  eStep <- list(xi=xi, M=matrix(0, nrow(data$Y), q), S=matrix(1e-4, nrow(data$Y), q))
  return(list(mStep=mStep, eStep=eStep, reg=reg, pca=pca))
}
OracleZiPLN <- function(data, latent){
  pca <- prcomp(latent$Z, rank=ncol(latent$W))
  mStep <- list(gamma=as.vector(glm(as.vector(latent$U) ~ -1 + data$X, family='binomial')$coef), 
                beta=as.vector(glm(as.vector(latent$Yall) ~ -1 + data$X + offset(as.vector(latent$Z)), 
                                   family='poisson')$coef), 
                C=pca$rotation%*%diag(pca$sdev[1:ncol(latent$W)]))
  eStep <- NULL
  return(list(mStep=mStep, eStep=eStep, pca=pca))
}

################################################################################
# ELBO
ELBOi <- function(datai, mStep, eStepi){
  nui <- as.vector(datai$Xi%*%mStep$gamma)
  mui <- as.vector(datai$Xi%*%mStep$beta)
  Ai <- exp(mui + as.vector(mStep$C%*%eStepi$mi) + 0.5*diag(mStep$C%*%diag(eStepi$Si)%*%t(mStep$C)))
  elboi <- sum(nui*eStepi$xii - log(1 + exp(nui)))
  elboi <- elboi - 0.5*(sum(eStepi$mi^2) + sum(eStepi$Si))
  elboi <- elboi + sum(datai$Omegai*eStepi$xii * 
                         (-Ai + datai$Yi*(mui + mStep$C%*%eStepi$mi) - datai$logFactYi))
  # Excludes xii=0 or 1, before computing the entropy
  xii <- eStepi$xii[which(eStepi$xii*(1-eStepi$xii) > 0)]
  elboi <- elboi + sum(xii * log(xii/(1 - xii)) + log(1 - xii))
  elboi <- elboi + 0.5*(length(eStepi$mi) + sum(log(eStepi$Si)))
  return(elboi)
}
ELBO <- function(data, mStep, eStep){
  sum(sapply(1:nrow(data$Y), function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], 
                  Omegai=data$Omega[i, ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=eStep$S[i, ])
    ELBOi(datai=datai, mStep=mStep, eStepi=eStepi)
  }))
}

################################################################################
# VE step
ElboMi <- function(mi, datai, mStep, eStepi){
  eStepi$mi <- mi
  ELBOi(datai=datai, mStep=mStep, eStepi=eStepi)
}
ElboGradMi <- function(mi, datai, mStep, eStepi){
  mui <- as.vector(datai$Xi%*%mStep$beta)
  Ai <- exp(mui + mi%*%t(mStep$C) + 0.5*diag(mStep$C%*%diag(eStepi$Si)%*%t(mStep$C)))
  as.vector(-mi + (datai$Omegai*eStepi$xii*(datai$Yi - Ai))%*%mStep$C)
}
ElboSi <- function(Si, datai, mStep, eStepi){
  eStepi$Si <- Si
  ELBOi(datai=datai, mStep=mStep, eStepi=eStepi)
}
ElboGradSi <- function(Si, datai, mStep, eStepi){
  mui <- as.vector(datai$Xi%*%mStep$beta)
  Ai <- exp(mui + eStepi$mi%*%t(mStep$C) + 0.5*diag(mStep$C%*%diag(Si)%*%t(mStep$C)))
  as.vector(0.5*(1/Si - 1 - (datai$Omegai*eStepi$xii*Ai)%*%(mStep$C^2)))
}
ElboMiSi <- function(miSi, datai, mStep, eStepi){
  q <- length(eStepi$mi)
  eStepi$mi <- miSi[1:q]; eStepi$Si <- miSi[q+(1:q)]
  ELBOi(datai=datai, mStep=mStep, eStepi=eStepi)
}
ElboGradMiSi <- function(miSi, datai, mStep, eStepi){
  q <- length(eStepi$mi)
  mi <- eStepi$mi <- miSi[1:q]; Si <- eStepi$Si <- miSi[q+(1:q)]
  c(ElboGradMi(mi, datai, mStep, eStepi), ElboGradSi(Si, datai, mStep, eStepi))
}


VEstep <- function(data, mStep, eStep, tolXi=1e-4, tolS=1e-4){
  n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X); q <- ncol(eStep$M)
  nuMuA <- NuMuA(data=data, mStep=mStep, eStep=eStep)
  nu <- nuMuA$nu; mu <- nuMuA$mu; A <- nuMuA$A
  xi <- plogis(nu - data$Omega*A)
  xi <- (xi + tolXi) / (1 + 2*tolXi)
  xi[which(data$Omega*data$Y>0)] <- 1
  eStepTmp <- list(xi=xi, M=eStep$M, S=eStep$S)
  if(ELBO(data=data, mStep=mStep, eStep=eStepTmp) > ELBO(data=data, mStep=mStep, eStep=eStep)){
    eStep$xi <- xi
  }
  # M <- t(sapply(1:n, function(i){
  #   datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], 
  #                 Omegai=data$Omega[i, ], logFactYi=data$logFactY[i, ])
  #   eStepi <- list(xii=eStep$xi[i, ], mi=NULL, Si=eStep$S[i, ])
  #   optim(par=eStep$M[i, ], fn=ElboMi, gr=ElboGradMi, 
  #         datai=datai, mStep=mStep, eStepi=eStepi, 
  #         method='BFGS', control=list(fnscale=-1))$par
  # }))
  # eStepTmp <- list(xi=eStep$xi, M=M, S=eStep$S)
  # if(ELBO(data=data, mStep=mStep, eStep=eStepTmp) > ELBO(data=data, mStep=mStep, eStep=eStep)){
  #   eStep$M <- M
  # }
  # S <- t(sapply(1:n, function(i){
  #   datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], 
  #                 Omegai=data$Omega[i, ], logFactYi=data$logFactY[i, ])
  #   eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=eStep$S[i, ])
  #   optim(par=eStepi$Si, fn=ElboSi, gr=ElboGradSi, 
  #         datai=datai, mStep=mStep, eStepi=eStepi, 
  #         method='L-BFGS-B', control=list(fnscale=-1), lower=rep(tolS, q))$par
  # }))
  # eStepTmp <- list(xi=eStep$xi, M=eStep$M, S=S)
  # if(ELBO(data=data, mStep=mStep, eStep=eStepTmp) > ELBO(data=data, mStep=mStep, eStep=eStep)){
  #   eStep$S <- S
  # }
  MS <- t(sapply(1:n, function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], 
                  Omegai=data$Omega[i, ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=eStep$S[i, ])
    optim(par=c(eStep$M[i, ], eStep$S[i, ]), fn=ElboMiSi, gr=ElboGradMiSi, 
          datai=datai, mStep=mStep, eStepi=eStepi, 
          method='L-BFGS-B', control=list(fnscale=-1), lower=c(rep(-Inf, q), rep(tolS, q)))$par
  }))
  eStepTmp <- list(xi=eStep$xi, M=MS[, 1:q], S=MS[, q+(1:q)])
  if(ELBO(data=data, mStep=mStep, eStep=eStepTmp) > ELBO(data=data, mStep=mStep, eStep=eStep)){
    eStep$M <- MS[, 1:q]; eStep$S <- MS[, q+(1:q)]
  }
  return(eStep)
}

################################################################################
# M step
ElboGamma <- function(gamma, data, mStep, eStep){
  mStep$gamma <- gamma
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboGradGamma <- function(gamma, data, mStep, eStep){
  nu <- matrix(data$X%*%gamma, nrow(data$Y), ncol(data$Y))
  probU <- plogis(nu)
  as.vector(t(data$X)%*%as.vector(eStep$xi - probU))
}
ElboBeta <- function(beta, data, mStep, eStep){
  mStep$beta <- beta
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboGradBeta <- function(beta, data, mStep, eStep){
  as.vector(t(data$X) %*% as.vector(data$Omega * eStep$xi * 
                                      (data$Y - NuMuA(data=data, mStep=mStep, eStep=eStep)$A)))
}
ElboC <- function(vecC, data, mStep, eStep){
  mStep$C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboGradC <- function(vecC, data, mStep, eStep){
  mStep$C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
  A <- NuMuA(data=data, mStep=mStep, eStep=eStep)$A
  as.vector(t(data$Omega*eStep$xi*(data$Y - A))%*%eStep$M - (t(data$Omega*eStep$xi*A)%*%eStep$S)*mStep$C)
}
Mstep <- function(data, mStep, eStep){
  beta <- optim(par=mStep$beta, fn=ElboBeta, gr=ElboGradBeta, data=data, mStep=mStep, eStep=eStep, 
                method='BFGS', control=list(fnscale=-1))$par
  mStep$beta <- beta
  mu <- matrix(data$X%*%mStep$beta, nrow(data$Y), ncol(data$Y))
  vecC <- optim(par=as.vector(mStep$C), fn=ElboC, gr=ElboGradC, data=data, mStep=mStep, eStep=eStep, 
                method='BFGS', control=list(fnscale=-1))$par
  mStep$C <- C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
  gamma <- optim(par=mStep$gamma, fn=ElboGamma, gr=ElboGradGamma, data=data, mStep=mStep, eStep=eStep, 
                 method='BFGS', control=list(fnscale=-1))$par
  mStep$gamma <- gamma
  return(mStep)
}

################################################################################
# VEM
VemZiPLN <- function(data, init, tol=1e-4, iterMax=1e3, tolXi=1e-4, tolS=1e-4, plot=TRUE){
  # init <- InitZiPLN(data); tol=1e-4; iterMax=1e3; tolXi=1e-5; tolS=1e-5
  mStep <- init$mStep; eStep <- init$eStep
  elboPath <- rep(NA, iterMax)
  diff <- 2*tol; iter <- 1
  elboPath[iter] <- ELBO(data=data, mStep=mStep, eStep=eStep)
  cat(iter, ':', elboPath[iter], diff)
  while((iter < iterMax) & (diff > tol)){
    iter <- iter+1
    mStepNew <- Mstep(data=data, mStep=mStep, eStep=eStep)
    eStepNew <- VEstep(data=data, mStep=mStepNew, eStep=eStep, tolXi=tolXi, tolS=tolS)
    elboPath[iter] <- ELBO(data=data, mStep=mStepNew, eStep=eStepNew)
    if(elboPath[iter] < elboPath[iter-1]){
      cat('', iter, ':', 'E-1', ELBO(data=data, mStep=mStep, eStep=eStep), 
          'M', ELBO(data=data, mStep=mStepNew, eStep=eStep), 
          'E', ELBO(data=data, mStep=mStepNew, eStep=eStepNew), '\n')
    }
    diff <- max(max(abs(eStepNew$M - eStep$M)), max(abs(mStepNew$C - mStep$C)))
    mStep <- mStepNew; eStep <- eStepNew
    if(iter%%round(sqrt(iterMax))==0){
      if(plot){plot(elboPath[1:iter], type='b', xlab='iter', ylab='elbo', 
                    ylim=quantile(elboPath[1:iter], probs=c(0.1, 1), na.rm=TRUE))}
      cat(' /', iter, ':', elboPath[iter], diff)
    }#else{cat('', iter)}
  }
  cat(' /', iter, ':', elboPath[iter], diff, '\n')
  pred <- NuMuA(data=data, mStep=mStep, eStep=eStep)
  pred$Yhat <- eStep$xi * pred$A
  elboPath <- elboPath[1:iter]
  return(list(mStep=mStep, eStep=eStep, pred=pred, iter=iter, elboPath=elboPath, elbo=elboPath[iter]))
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
  condPi <- margPi <- plogis(nuMuA$nu)
  margLambda <- exp(nuMuA$nu + rep(1, n)%o%diag(sigma)/2)
  condLambda <- nuMuA$A
  nuMuA <-NuMuA(data=data, mStep=fit$mStep, eStep=fit$eStep)
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
