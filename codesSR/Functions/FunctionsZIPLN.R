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
  reg <- lm(as.vector(log(1+data$Y)) ~ -1 + data$X)
  pca <- prcomp(matrix(reg$residuals, n, p), rank=q)
  mStep <- list(gamma=as.vector(glm(as.vector(1*(data$Y > 0)) ~ -1 + data$X, family='binomial')$coef),
                beta=reg$coef,
                C=pca$rotation %*% diag(pca$sdev[1:q]))
  eStep <- list(xi=matrix(mean(data$Y > 0), n, p), M=matrix(0, n, q), S=matrix(1e-4, n, q))
  return(list(mStep=mStep, eStep=eStep, reg=reg, pca=pca))
}
InitZiPLN <- function(data){
  reg <- lm(as.vector(log(1+data$Y)) ~ -1 + data$X)
  pca <- prcomp(matrix(reg$residuals, n, p), rank=q)
  zip <- EmZIP(X=data$X, Y=as.vector(data$Y))
  mStep <- list(gamma=zip$gamma, beta=zip$beta, 
                C=pca$rotation %*% diag(pca$sdev[1:q]))
  eStep <- list(xi=matrix(mean(data$Y > 0), n, p), M=matrix(0, n, q), S=matrix(1e-4, n, q))
  eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep)
  return(list(mStep=mStep, eStep=eStep, reg=reg, pca=pca))
}
OracleZiPLN <- function(sim){
  pca <- prcomp(sim$Z, rank=ncol(sim$W))
  mStep <- list(gamma=as.vector(glm(as.vector(sim$U) ~ -1 + sim$X, family='binomial')$coef), 
                beta=as.vector(glm(as.vector(sim$Yall) ~ -1 + sim$X + offset(as.vector(sim$Z)), 
                                   family='poisson')$coef), 
                C=pca$rotation%*%diag(pca$sdev[1:ncol(sim$W)]))
  eStep <- NULL
  return(list(mStep=mStep, eStep=eStep, pca=pca))
}

################################################################################
# ELBO
ELBOi <- function(datai, mStep, eStepi){
  nui <- as.vector(datai$Xi%*%mStep$gamma)
  mui <- as.vector(datai$Xi%*%mStep$beta)
  Ai <- exp(mui + as.vector(mStep$C%*%eStepi$mi) + 0.5*diag(mStep$C%*%diag(eStepi$Si)%*%t(mStep$C)))
  elboi <- sum(nui*eStepi$xii - log(1 + exp(nui))) # 0K 
  elboi <- elboi - 0.5*(sum(eStepi$mi^2) + sum(eStepi$Si)) #OK 
  elboi <- elboi + sum(eStepi$xii * (-Ai + datai$Yi*(mui + mStep$C%*%eStepi$mi) - datai$logFactYi)) # OK 
  # Excludes xii=0 or 1, before computing the entropy
  xii <- eStepi$xii[which(eStepi$xii*(1-eStepi$xii) > 0)] # OK 
  elboi <- elboi + sum(xii * log(xii/(1 - xii)) + log(1 - xii)) # OK 
  elboi <- elboi + 0.5*(length(eStepi$mi) + sum(log(eStepi$Si))) # OK
  return(elboi)
}
ELBO <- function(data, mStep, eStep){
  sum(sapply(1:nrow(data$Y), function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=eStep$S[i, ])
    ELBOi(datai=datai, mStep=mStep, eStepi=eStepi)
  }))
}


# ELBO
ELBOi <- function(datai, mStep, eStepi){
  nui <- as.vector(datai$Xi%*%mStep$gamma)
  mui <- as.vector(datai$Xi%*%mStep$beta)
  Ai <- exp(mui + as.vector(mStep$C%*%eStepi$mi) + 0.5*diag(mStep$C%*%diag(eStepi$Si)%*%t(mStep$C)))
  elboi <- sum(nui*eStepi$xii - log(1 + exp(nui))) # 0K 
  elboi <- elboi - 0.5*(sum(eStepi$mi^2) + sum(eStepi$Si)) #OK 
  elboi <- elboi + sum(eStepi$xii * (-Ai + datai$Yi*(mui + mStep$C%*%eStepi$mi) - datai$logFactYi)) # OK 
  # Excludes xii=0 or 1, before computing the entropy
  xii <- eStepi$xii[which(eStepi$xii*(1-eStepi$xii) > 0)] # OK 
  elboi <- elboi + sum(xii * log(xii/(1 - xii)) + log(1 - xii)) # OK 
  elboi <- elboi + 0.5*(length(eStepi$mi) + sum(log(eStepi$Si))) # OK
  return(elboi)
}
ELBO <- function(data, mStep, eStep){
  sum(sapply(1:nrow(data$Y), function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=eStep$S[i, ])
    ELBOi(datai=datai, mStep=mStep, eStepi=eStepi)
  }))
}

ELBOSophiei <- function(datai, mStep, eStepi){
  nui <- as.vector(datai$Xi%*%mStep$gamma)
  mui <- as.vector(datai$Xi%*%mStep$beta)
  Zi <- mui + as.vector(mStep$C%*%eStepi$mi)
  q <- dim(mStep$C)[2]
  p <- dim(mStep$C)[1]
  Ai <- exp(Zi + 0.5*diag(mStep$C%*%diag(eStepi$Si,q,q)%*%t(mStep$C)))
  S1i <- sum(nui*eStepi$xii - log(1 + exp(nui)))
  S2i <- - 0.5*(sum(eStepi$mi^2) + sum(eStepi$Si))
  S3i <- sum(eStepi$xii * (-Ai + datai$Yi*(mui + mStep$C%*%eStepi$mi) - datai$logFactYi))
  xii <- eStepi$xii[which(eStepi$xii*(1-eStepi$xii) > 0)]
  S4i <- sum(xii * log(xii/(1 - xii)) + log(1 - xii)) 
  S5i <- 0.5*(length(eStepi$mi) + sum(log(eStepi$Si)))
  elboi = S1i + S2i + S3i + S4i+ S5i
  return(c(S1i,S2i,S3i,S4i,S5i,elboi))
}
ELBOSophie <- function(data, mStep, eStep){
  US <- sapply(1:nrow(data$Y), function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=eStep$S[i, ])
    ELBOSophiei(datai=datai, mStep=mStep, eStepi=eStepi)
  })
  return(apply(US,1,sum))
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
  as.vector(-mi + (eStepi$xii*(datai$Yi - Ai))%*%mStep$C)
}
ElboSi <- function(Si, datai, mStep, eStepi){
  eStepi$Si <- Si
  ELBOi(datai=datai, mStep=mStep, eStepi=eStepi)
}
ElboGradSi <- function(Si, datai, mStep, eStepi){
  mui <- as.vector(datai$Xi%*%mStep$beta)
  Ai <- exp(mui + eStepi$mi%*%t(mStep$C) + 0.5*diag(mStep$C%*%diag(Si)%*%t(mStep$C)))
  as.vector(0.5*(1/Si - 1 - (eStepi$xii*Ai)%*%(mStep$C^2)))
}
# ElboM <- function(Mvec, data, mStep, eStep){
#   eStep$M <- matrix(Mvec, nrow(data$Y), ncol(eStep$M))
#   ELBO(data=data, mStep=mStep, eStep=eStep)
# }
# ElboGradM <- function(Mvec, data, mStep, eStep){
#   M <- matrix(Mvec, nrow(data$Y), ncol(eStep$M))
#   as.vector(t(sapply(1:nrow(data$Y), function(i){
#     datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
#     eStepi <- list(xii=eStep$xi[i, ], mi=NULL, Si=eStep$S[i, ])
#     ElboGradMi(mi=M[i, ], datai=datai, mStep=mStep, eStepi=eStepi)
#   })))
# }
# ElboS <- function(Svec, datai, mStep, eStepi){
#   eStep$S <- matrix(Svec, nrow(data$Y), ncol(eStep$M))
#   ELBO(data=data, mStep=mStep, eStep=eStep)
# }
# ElboGradS <- function(S, data, mStep, eStep){
#   S <- matrix(Svec, nrow(data$Y), ncol(eStep$M))
#   as.vector(t(sapply(1:nrow(data$Y), function(i){
#     datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
#     eStepi <- list(xii=eStep$xi[i, ], mi=NULL, Si=eStep$S[i, ])
#     ElboGradSi(Si=S[i, ], datai=datai, mStep=mStep, eStepi=eStepi)
#   })))
# }

VEstep <- function(data, mStep, eStep, tolXi=1e-4, tolS=1e-4){
  n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X); q <- ncol(eStep$M)
  nuMuA <- NuMuA(data=data, mStep=mStep, eStep=eStep)
  nu <- nuMuA$nu; mu <- nuMuA$mu; A <- nuMuA$A
  xi <- plogis(nu - A)
  xi <- (xi + tolXi) / (1 + 2*tolXi)
  xi[which(data$Y>0)] <- 1
  eStepTmp <- list(xi=xi, M=eStep$M, S=eStep$S)
  if(ELBO(data=data, mStep=mStep, eStep=eStepTmp) > ELBO(data=data, mStep=mStep, eStep=eStep)){
    eStep$xi <- xi
  }
  M <- t(sapply(1:n, function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=NULL, Si=eStep$S[i, ])
    optim(par=eStep$M[i, ], fn=ElboMi, gr=ElboGradMi, 
          datai=datai, mStep=mStep, eStepi=eStepi, 
          method='BFGS', control=list(fnscale=-1))$par
  }))
  eStepTmp <- list(xi=eStep$xi, M=M, S=eStep$S)
  if(ELBO(data=data, mStep=mStep, eStep=eStepTmp) > ELBO(data=data, mStep=mStep, eStep=eStep)){
    eStep$M <- M
  }
  S <- t(sapply(1:n, function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=eStep$S[i, ])
    optim(par=eStepi$Si, fn=ElboSi, gr=ElboGradSi, 
          datai=datai, mStep=mStep, eStepi=eStepi, 
          method='L-BFGS-B', control=list(fnscale=-1), lower=rep(tolS, q))$par
  }))
  eStepTmp <- list(xi=eStep$xi, M=eStep$M, S=S)
  if(ELBO(data=data, mStep=mStep, eStep=eStepTmp) > ELBO(data=data, mStep=mStep, eStep=eStep)){
    eStep$S <- S
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
  as.vector(t(data$X) %*% as.vector(eStep$xi*(data$Y - NuMuA(data=data, mStep=mStep, eStep=eStep)$A)))
}
ElboC <- function(vecC, data, mStep, eStep){
  mStep$C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboGradC <- function(vecC, data, mStep, eStep){
  mStep$C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
  A <- NuMuA(data=data, mStep=mStep, eStep=eStep)$A
  as.vector(t(eStep$xi*(data$Y - A))%*%eStep$M - (t(eStep$xi*A)%*%eStep$S)*mStep$C)
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
      plot(elboPath[1:iter], type='b', xlab='iter', ylim=quantile(elboPath[1:iter], probs=c(0.1, 1), na.rm=TRUE))
      cat(' /', iter, ':', elboPath[iter], diff)
    }#else{cat('', iter)}
  }
  cat('\n')
  pred <- NuMuA(data=data, mStep=mStep, eStep=eStep)
  pred$Yhat <- eStep$xi * pred$A
  elboPath <- elboPath[1:iter]
  return(list(mStep=mStep, eStep=eStep, pred=pred, iter=iter, elboPath=elboPath, elbo=elboPath[iter]))
}
