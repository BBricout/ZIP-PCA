# Functions ZI-PLN

################################################################################
# Reshape
Mstep2Theta <- function(mStep, n, d, p, q){
  as.vector(c(mStep$gamma, mStep$beta, as.vector(mStep$C)))
}
Theta2Mstep <- function(theta, n, d, p, q){
  list(gamma=theta[1:d], beta=theta[d+(1:d)], C=matrix(theta[2*d+(1:(p*q))], p, q))
}
Estep2Psi <- function(eStep, n, d, p, q){
  as.vector(c(as.vector(eStep$M), as.vector(eStep$S)))
}
Psi2Estep <- function(theta, n, d, p, q){
  list(M=matrix(psi[(1:(n*q))], n, q), S=matrix(psi[(n*q) + (1:(n*q))], n, q))
}
ComputeXi <- function(data, mStep, estep){
  nuMuA <- NuMuA(data=data, mStep=mStep, eStep=eStep)
  nu <- nuMuA$nu; mu <- nuMuA$mu; A <- nuMuA$A
  xi <- plogis(nu - data$Omega*A)
  xi <- (xi + tolXi) / (1 + 2*tolXi)
  xi[which(data$Omega*data$Y>0)] <- 1
  return(xi)
}

################################################################################
# ELBO
# mStep <- init$mStep; eStep <- init$eStep; 
# theta <- Mstep2Theta(mStep, n, d, p, q); psi <- Estep2Psi(eStep, n, d, p, q); 
# thetaPsi <- c(theta, psi)
ElboVecTheta <- function(theta, data, eStep){
  mStep <- Theta2Mstep(theta, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), q=ncol(eStep$M))
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboVecPsi <- function(psi, data, eStep){
  eStep <- Psi2Estep(psi, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), q=ncol(eStep$M))
  eStep$xi <- ComputeXi(data=data, mStep=mStep)
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboVecThetaPsi <- function(thetaPsi, data, n, d, p, q){
  theta <- thetaPsi[1:((2*d)+(p*q))]
  mStep <- Theta2Mstep(theta, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), q=ncol(eStep$M))
  psi <- thetaPsi[(2*d)+(p*q) + (1:(2*(n*q)))]
  eStep <- Psi2Estep(psi, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), q=ncol(eStep$M))
  eStep$xi <- ComputeXi(data=data, mStep=mStep)
  ELBO(data=data, mStep=mStep, eStep=eStep)
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
    eStepi <- list(xii=eStep$xi[i, ], mi=Mvec[(i-1)*q + (1:q)], Si=NULL)
    ElboGradSi(Si=Si, datai=datai, mStep=mStep, eStepi=eStepi)
  }))
}
ElboGradVecPsi <- function(psi, data, mStep){
  eStep <- Psi2Estep(psi, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), q=ncol(eStep$M))
  eStep$xi <- ComputeXi(data=data, mStep=mStep)
  c(ElboGradVecM(Mvec=as.vector(t(eStep$M)), data=data, mStep=mStep, eStep=eStep), 
    ElboGradVecS(Svec=as.vector(t(eStep$S)), data=data, mStep=mStep, eStep=eStep))
}
ElboGradVecThetaPsi <- function(thetaPsi, data, q){
  n <- nrow(data$Y); d <- ncol(data$X); p <- ncol(data$Y)
  theta <- thetaPsi[1:((2*d)+(p*q))]
  mStep <- Theta2Mstep(theta, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), q=ncol(eStep$M))
  psi <- thetaPsi[(2*d)+(p*q) + (1:(2*(n*q)))]
  eStep <- Psi2Estep(psi, n=nrow(data$Y), d=ncol(data$X), p=ncol(data$Y), q=ncol(eStep$M))
  eStep$xi <- ComputeXi(data=data, mStep=mStep)
  c(ElboGradVecTheta(theta, data, eStep), ElboGradVecPsi(psi, data, mStep))
}

  