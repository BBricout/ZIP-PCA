# Functions ZI-PLN

################################################################################
# ELBO
# mStep <- init$mStep; eStep <- init$eStep; 
# theta <- Mstep2Theta(mStep, n, d, p, q); psi <- Estep2Psi(eStep, n, d, p, q); 
# thetaPsi <- c(theta, psi)
ElboVecTheta <- function(theta, data, eStep){
  mStep <- Theta2Mstep(theta, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), q=ncol(eStep$M))
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboVecPsi <- function(psi, data, mStep, tolXi=1e-4){
  eStep <- Psi2Estep(psi, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), 
                     q=round(length(psi)/2/nrow(data$Y)))
  eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=tolXi)
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboVecThetaPsi <- function(thetaPsi, data, q, tolXi=1e-4){
  n <- nrow(data$Y); d <- ncol(data$X); p <- ncol(data$Y)
  theta <- thetaPsi[1:((2*d)+(p*q))]
  mStep <- Theta2Mstep(theta, n=n, d=d, p=p, q=q)
  psi <- thetaPsi[(2*d)+(p*q) + (1:(2*(n*q)))]
  eStep <- Psi2Estep(psi, n=n, d=d, p=p, q=q)
  eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=tolXi)
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
NegElboVecThetaPsi <- function(thetaPsi, data, q, tolXi=1e-4){
  -ElboVecThetaPsi(thetaPsi=thetaPsi, data=data, q=q, tolXi=tolXi)
}

################################################################################
# ELBO grad
# mStep <- init$mStep; eStep <- init$eStep; 
# theta <- Mstep2Theta(mStep, n, d, p, q); psi <- Estep2Psi(eStep, n, d, p, q); 
# thetaPsi <- c(theta, psi)
ElboGradVecTheta <- function(theta, data, eStep){
  mStep <- Theta2Mstep(theta, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), q=ncol(eStep$M))
  c(ElboGradGamma(mStep$gamma, data=data, mStep=mStep, eStep=eStep), 
    ElboGradBeta(mStep$beta, data=data, mStep=mStep, eStep=eStep), 
    ElboGradC(as.vector(mStep$C), data=data, mStep=mStep, eStep=eStep))
}
ElboGradVecM <- function(Mvec, data, mStep, eStep){
  n <- nrow(data$Y); q <- ncol(mStep$C)
  as.vector(sapply(1:n, function(i){
    mi <- Mvec[(i-1)*q + (1:q)]
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ],
                  Omegai=data$Omega[i, ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=NULL, Si=eStep$S[i, ])
    ElboGradMi(mi=mi, datai=datai, mStep=mStep, eStepi=eStepi)
    }))
}
ElboGradVecS <- function(Svec, data, mStep, eStep){
  n <- nrow(data$Y); q <- ncol(mStep$C)
  as.vector(sapply(1:n, function(i){
    Si <- Svec[(i-1)*q + (1:q)]
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ],
                  Omegai=data$Omega[i, ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=NULL)
    ElboGradSi(Si=Si, datai=datai, mStep=mStep, eStepi=eStepi)
  }))
}
ElboGradVecPsi <- function(psi, data, mStep, tolXi=1e-4){
  eStep <- Psi2Estep(psi, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), 
                     q=round(length(psi)/2/nrow(data$Y)))
  eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=tolXi)
  c(ElboGradVecM(Mvec=as.vector(t(eStep$M)), data=data, mStep=mStep, eStep=eStep), 
    ElboGradVecS(Svec=as.vector(t(eStep$S)), data=data, mStep=mStep, eStep=eStep))
}
ElboGradVecThetaPsi <- function(thetaPsi, data, q, tolXi=1e-4){
  n <- nrow(data$Y); d <- ncol(data$X); p <- ncol(data$Y)
  theta <- thetaPsi[1:((2*d)+(p*q))]
  mStep <- Theta2Mstep(theta, n=n, d=d, p=p, q=q)
  psi <- thetaPsi[(2*d)+(p*q) + (1:(2*(n*q)))]
  eStep <- Psi2Estep(psi, n=n, d=d, p=p, q=q)
  eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=tolXi)
  c(ElboGradVecTheta(theta, data, eStep), ElboGradVecPsi(psi, data, mStep))
}
NegElboGradVecThetaPsi <- function(thetaPsi, data, q, tolXi=1e-4){
  -ElboGradVecThetaPsi(thetaPsi=thetaPsi, data=data, q=q, tolXi=tolXi)
}

################################################################################
# ELBO + Grad : log(S) versions
ElboLogSVecThetaPsi <- function(thetaPsiLogS, data, q, tolXi=1e-4){
  n <- nrow(data$Y); d <- ncol(data$X); p <- ncol(data$Y)
  theta <- thetaPsiLogS[1:((2*d)+(p*q))]
  mStep <- Theta2Mstep(theta, n=n, d=d, p=p, q=q)
  psiLogS <- thetaPsiLogS[(2*d)+(p*q) + (1:(2*(n*q)))]
  eStep <- Psi2Estep(psiLogS, n=n, d=d, p=p, q=q)
  eStep$S <- exp(eStep$S)
  eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=tolXi)
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
NegElboLogSVecThetaPsi <- function(thetaPsiLogS, data, q, tolXi=1e-4){
  -ElboLogSVecThetaPsi(thetaPsiLogS=thetaPsiLogS, data=data, q=q, tolXi=tolXi)
}

ElboLogSGradVecThetaPsi <- function(thetaPsiLogS, data, q, tolXi=1e-4){
  n <- nrow(data$Y); d <- ncol(data$X); p <- ncol(data$Y)
  theta <- thetaPsiLogS[1:((2*d)+(p*q))]
  mStep <- Theta2Mstep(theta, n=n, d=d, p=p, q=q)
  psi <- psiLogS <- thetaPsiLogS[(2*d)+(p*q) + (1:(2*(n*q)))]
  psi[(n*q)+(1:(n*q))] <- exp(psi[(n*q)+(1:(n*q))])
  eStep <- Psi2Estep(psi, n=n, d=d, p=p, q=q)
  eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=tolXi)
  gradTheta <- ElboGradVecTheta(theta, data, eStep)
  gradPsi <- ElboGradVecPsi(psi, data, mStep) * c(rep(1, n*q), psi[(n*q)+(1:(n*q))])
  return(c(gradTheta, gradPsi))
}
NegElboLogSGradVecThetaPsi <- function(thetaPsiLogS, data, q, tolXi=1e-4){
  -ElboLogSGradVecThetaPsi(thetaPsiLogS=thetaPsiLogS, data=data, q=q, tolXi=tolXi)
}

################################################################################
# VEM
VemZiPLNvec <- function(data, init, tol=1e-5, iterMax=5e3, iterMaxOptim=1e1, tolXi=1e-4, tolS=1e-4, plot=TRUE){
  # tol=1e-4; iterMax=1e3; tolXi=1e-4; tolS=1e-4; plot=TRUE
  n <- nrow(data$Y); d <- ncol(data$X); q <- ncol(init$eStep$M); p <- ncol(data$Y)
  mStep <- init$mStep; theta <- Mstep2Theta(mStep, n=n, d=d, p=p, q=q)
  eStep <- init$eStep; psi <- Estep2Psi(eStep, n=n, d=d, p=p, q=q)
  elboPath <- rep(NA, iterMax)
  iter <- 1; diff <- 2*tol
  elboPath[iter] <- ElboVecThetaPsi(thetaPsi=c(theta, psi), data=data, q=q)
  while((iter < iterMax) & (diff > tol)){
    iter <- iter+1
    # # VE + M
    # fit <- optim(par=c(theta, psi), fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=data,
    #              q=q, method='L-BFGS-B', control=list(fnscale=-1),
    #              lower=c(rep(-Inf, ((2*d) + (p*q) + (n*q))), rep(tolS, (n*q))))
    # thetaNew <- fit$par[1:((2*d)+(p*q))]; psiNew <- fit$par[-(1:((2*d)+(p*q)))]
    # VE
    fitVE <- optim(par=psi, fn=ElboVecPsi, gr=ElboGradVecPsi, data=data, mStep=mStep,
                   method='L-BFGS-B', control=list(fnscale=-1, maxit=iterMaxOptim),
                   lower=c(rep(-Inf, (n*q)), rep(tolS, (n*q))))
    psiNew <- fitVE$par; eStepNew <- Psi2Estep(psiNew, n=n, d=d, p=p, q=q)
    eStepNew$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStepNew, tolXi=tolXi)
    # M
    fitM <- optim(par=theta, fn=ElboVecTheta, gr=ElboGradVecTheta, data=data, eStep=eStepNew, 
                  method='BFGS', control=list(fnscale=-1, maxit=iterMaxOptim))
    thetaNew <- fitM$par; mStepNew <- Theta2Mstep(thetaNew, n=n, d=d, p=p, q=q)
    #   
    elboPath[iter] <- ElboVecThetaPsi(thetaPsi=c(thetaNew, psiNew), data=data, q=q)
    diff <- max(abs(c(thetaNew, psiNew) - c(theta, psi)))  
    theta <- thetaNew; mStep <- Theta2Mstep(theta, n=n, d=d, p=p, q=q)
    psi <- psiNew; eStep <- Psi2Estep(psi, n=n, d=d, p=p, q=q)
    if(iter%%round(sqrt(iterMax))==0){
      cat(q, iter, ':', diff, fitVE$value, fitM$value, elboPath[iter], '\n')
      plot(elboPath[1:iter], type='b', ylim=quantile(elboPath[1:iter], prob=c(0.1, 1)))
    }
  }
  elboPath <- elboPath[1:iter]
  pred <- NuMuA(data=data, mStep=mStep, eStep=eStep)
  pred$Yhat <- eStep$xi * pred$A
  return(list(mStep=mStep, eStep=eStep, pred=pred, iter=iter, elboPath=elboPath, elbo=elboPath[iter]))
}

