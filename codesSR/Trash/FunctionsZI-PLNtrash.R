ElboM <- function(vecM, data, mStep, eStep){
  eStep$M <- matrix(vecM, nrow(data$Y), ncol(mStep$C))
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboS <- function(vecS, data, mStep, eStep){
  eStep$S <- matrix(vecS, nrow(data$Y), ncol(mStep$C))
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboMS <- function(vecMS, data, mStep, eStep){
  n <- nrow(data$Y); q <- ncol(mStep$C)
  eStep$M <- matrix(vecMS[1:(p*q)], p, q)
  eStep$S <- matrix(vecMS[(p*q)+(1:(p*q))], p, q)
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboGradM <- function(vecM, data, mStep, eStep){
  M <- matrix(vecM, nrow(data$Y), ncol(mStep$C))
  sapply(1:n, function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=NULL, Si=eStep$S[i, ])
    ElboGradMi(mi=M[i, ], datai=datai, mStep=mStep, eStepi=eStepi)
    })
}
ElboGradS <- function(vecS, data, mStep, eStep){
  S <- matrix(vecS, nrow(data$Y), ncol(mStep$C))
  sapply(1:n, function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=NULL)
    ElboGradSi(Si=S[i, ], datai=datai, mStep=mStep, eStepi=eStepi)
  })
}
ElboGradMS <- function(vecMS, data, mStep, eStep){
  n <- nrow(data$Y); q <- ncol(mStep$C)
  vecM <- vecMS[1:(p*q)]; eStep$M <- matrix(vecM, p, q)
  vecS <- vecMS[(p*q)+(1:(p*q))]; eStep$S <- matrix(vecS, p, q)
  c(ElboGradM(vecM=vecM, data=data, mstep=mStep, eStep=eStep), 
    ElboGradS(vecS=vecS, data=data, mStep=mStep, eStep=eStep))
}

VEstepTrash <- function(data, mStep, eStep, tolXi=1e-4){
  n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X); q <- ncol(eStep$M)
  nu <- matrix(data$X%*%mStep$gamma, n, p)
  mu <- matrix(data$X%*%mStep$beta, n, p)
  A <- exp(mu + eStep$M%*%t(mStep$C) + 
             0.5*t(sapply(1:n, function(i){diag(mStep$C%*%diag(eStep$S[i, ])%*%t(mStep$C))})))
  xiLogit <- nu - A + data$Y*(eStep$M%*%t(mStep$C) + mu) - data$logFactY
  xi <- plogis(xiLogit)
  xi <- (xi + tolXi) / (1 + 2*tolXi)
  eStep$xi <- xi
  Mvec <- optim(par=as.vector(eStep$M), fn=ElboM, gr=ElboGradM, 
          data=data, mStep=mStep, eStep=eStep, 
          method='BFGS', control=list(fnscale=-1))$par
  M <- eStep$M <- matrix(Mvec, n, q)
  Svec <- optim(par=as.vector(eStep$S), fn=ElboS, gr=ElboGradS, 
                data=data, mStep=mStep, eStep=eStep, 
                method='L-BFGS-B', control=list(fnscale=-1), lower=rep(1e-4, n*q))$par
  S <- matrix(Svec, n, q)
  return(list(xi=xi, M=M, S=S))
}

################################################################################
# M step
ElboGammaBetaC <- function(gammaBetaC, data, mStep, eStep){
  n <- nrow(data$Y); p <- ncol(data$Y); q <- ncol(eStep$M)
  mStep$gamma <- gammaBetaC[1:d]
  mStep$beta <- gammaBetaC[d+(1:d)]
  vecC <- gammaBetaC[2*d+(1:p*q)]
  mStep$C <- matrix(vecC, p, q)
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboGradGammaBetaC <- function(gammaBetaC, data, mStep, eStep){
  n <- nrow(data$Y); p <- ncol(data$Y); q <- ncol(eStep$M)
  mStep$gamma <- gammaBetaC[1:d]
  mStep$beta <- gammaBetaC[d+(1:d)]
  vecC <- gammaBetaC[2*d+(1:p*q)]
  mStep$C <- matrix(vecC, p, q)
  c(ElboGradGamma(gamma=mStep$gamma, data=data, mStep=mStep, eStep=eStep), 
    ElboGradBeta(beta=mStep$beta, data=data, mStep=mStep, eStep=eStep), 
    ElboGradC(vecC=vecC, data=data, mStep=mStep, eStep=eStep))
}

MstepTrash <- function(data, mStep, eStep){
  n <- nrow(data$Y); p <- ncol(data$Y); q <- ncol(eStep$M)
  gammaBetaC <- optim(par=c(mStep$gamma, mStep$beta, as.vector(mStep$C)), 
                      fn=ElboGammaBetaC, gr=ElboGradGammaBetaC, 
                      data=data, mStep=mStep, eStep=eStep, 
                      method='BFGS', control=list(fnscale=-1))$par
  return(list(gamma=gammaBetaC[1:d], 
              beta=gammaBetaC[d+(1:d)], 
              C=matrix(gammaBetaC[2*d+(1:p*q)], p, q)))
}

MstepC <- function(data, mStep, eStep){
  vecC <- optim(par=as.vector(mStep$C), fn=ElboC, gr=ElboGradC, data=data, mStep=mStep, eStep=eStep,
                method='BFGS', control=list(fnscale=-1))$par
  mStep$C <- C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
  # svdC <- svd(matrix(vecC, ncol(data$Y), ncol(eStep$M)))
  # mStep$C <- C <- svdC$u%*%diag(svdC$d)
  return(list(gamma=mStep$gamma, beta=mStep$beta, C=C))
}
