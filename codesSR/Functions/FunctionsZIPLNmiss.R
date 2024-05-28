# Functions ZI-PLN

################################################################################
# Simul
NuMuA <- function(data, mStep, eStep){
  nu <- matrix(data$X%*%mStep$gamma, n, p)
  mu <- matrix(data$X%*%mStep$beta, n, p)
  A <- exp(mu + eStep$M%*%t(mStep$C) + 
             0.5*t(sapply(1:n, function(i){diag(mStep$C%*%diag(eStep$S[i, ])%*%t(mStep$C))})))
  return(list(nu=nu, mu=mu, A=A))
}

################################################################################
# Init
InitZiPLNold <- function(data){
  pres <- which(data$Omega==1)
  reg <- lm(as.vector(log(1+data$Y)[pres]) ~ -1 + data$X[pres, ])
  res <- matrix(0, nrow(data$Y), ncol(data$Y))
  res[pres] <- reg$residuals
  pca <- prcomp(res, rank=q)
  mStep <- list(gamma=as.vector(glm(as.vector(1*(data$Y > 0)) ~ -1 + data$X, family='binomial')$coef), 
                beta=reg$coef, 
                C=pca$rotation %*% diag(pca$sdev[1:q]))
  eStep <- list(xi=matrix(sum(data$Omega*(data$Y > 0))/sum(data$Omega), n, p), M=matrix(0, n, q), S=matrix(1e-4, n, q))
  return(list(mStep=mStep, eStep=eStep, reg=reg, pca=pca))
}
InitZiPLN <- function(data){
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
  xi <- matrix(plogis(data$X%*%mStep$gamma), n, p)
  xi[which(data$Y > 0)] <- 1
  eStep <- list(xi=xi, M=matrix(0, n, q), S=matrix(1e-4, n, q))
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
  mu <- matrix(data$X%*%mStep$beta, n, p)
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
VemZiPLN <- function(data, init, tol=1e-4, iterMax=1e3, tolXi=1e-4, tolS=1e-4){
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
      plot(elboPath[1:iter], type='b', xlab='iter', ylab='elbo', 
           ylim=quantile(elboPath[1:iter], probs=c(0.1, 1), na.rm=TRUE))
      cat(' /', iter, ':', elboPath[iter], diff)
    }#else{cat('', iter)}
  }
  cat('\n')
  pred <- NuMuA(data=data, mStep=mStep, eStep=eStep)
  pred$Yhat <- eStep$xi * pred$A
  elboPath <- elboPath[1:iter]
  return(list(mStep=mStep, eStep=eStep, pred=pred, iter=iter, elboPath=elboPath, elbo=elboPath[iter]))
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
  margCIl <- margCIu <- condCIl <- condCIu <- matrix(NA, n, p)
  for(i in 1:n){for(j in 1:p){
    margCIl[i, j] <- qzip(p=.025, psi=1-margPi[i, j], lambda=margPi[i, j])
    margCIu[i, j] <- qzip(p=.975, psi=1-margPi[i, j], lambda=margPi[i, j])
    condCIl[i, j] <- qzip(p=.025, psi=1-condPi[i, j], lambda=condPi[i, j])
    condCIu[i, j] <- qzip(p=.975, psi=1-condPi[i, j], lambda=condPi[i, j])
  }}
  return(list(nuMuA=nuMuA,
              marg=list(pi=margPi, lambda=margLambda,
                        esp=margEsp, var=margVar, lCI=margCIl, uCI=margCIu), 
              cond=list(pi=condPi, lambda=condLambda,
                        esp=condEsp, var=condVar, lCI=condCIl, uCI=condCIu)))
}
