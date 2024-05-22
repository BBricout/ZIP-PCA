# Functions ZI-PLN

################################################################################
# Simul
SimZiPLN <- function(n, p, d, q){
  X <- matrix(rnorm(n*p*d), n*p, d); X[, 1] <- 1
  ij <- cbind(rep(1:n, p), rep(1:p, each=n))
  # presence
  gamma <- rnorm(d)/sqrt(d); gamma[1] <- gamma[1] + 2
  nu <- matrix(X%*%gamma, n, p)
  probU <- plogis(nu)
  U <- matrix(rbinom(n*p, 1, probU), n, p)
  # latent
  W <- matrix(rnorm(n*q), n, q)
  C <- matrix(rnorm(p*q), p, q)
  Z <- W%*%t(C)
  # abundance
  beta <- rnorm(d)/sqrt(d); beta[1] <- beta[1] + 2
  mu <- matrix(X%*%beta, n, p)
  lambdaY <- exp(mu)
  Y <- matrix(rpois(n*p, lambdaY), n, p)
  return(list(X=X, Y=Y, ij=ij, U=U, W=W, Z=Z, gamma=gamma, beta=beta, gamma=gamma, C=C))
}

################################################################################
# Init
InitZiPLN <- function(data){
  reg <- lm(as.vector(log(1+data$Y)) ~ -1 + data$X)
  pca <- prcomp(matrix(reg$residuals, n, p), rank=q)
  pca <- pca$rotation %*% diag(pca$sdev[1:q])
  mStep <- list(gamma=as.vector(glm(as.vector(1*(data$Y > 0)) ~ -1 + data$X, family='binomial')$coef), 
                beta=as.vector(reg$coef), 
                C=pca)
  eStep <- list(xi=matrix(mean(data$Y > 0), n, p), M=matrix(0, n, q), S=true$S)
  return(list(mStep=mStep, eStep=eStep, reg=reg, pca=pca))
}
################################################################################
# ELBO
ELBOi <- function(datai, mStep, eStepi){
  nui <- as.vector(datai$Xi%*%mStep$gamma)
  mui <- as.vector(datai$Xi%*%mStep$beta)
  Ai <- exp(mui + as.vector(mStep$C%*%eStepi$mi) + 0.5*diag(mStep$C%*%diag(eStepi$Si)%*%t(mStep$C)))
  elboi <- sum(nui*eStepi$xii - log(1 + exp(nui)))
  elboi <- elboi - 0.5*(sum(eStepi$mi^2) + sum(eStepi$Si))
  elboi <- elboi + sum(eStepi$xii * (-Ai + datai$Yi*(mui + mStep$C%*%eStepi$mi) - datai$logFactYi))
  elboi <- elboi + sum(eStepi$xii * log(eStepi$xii/(1 - eStepi$xii)) + log(1 - eStepi$xii))
  elboi <- elboi + 0.5*(length(eStepi$mi) + sum(log(eStepi$Si)))
  return(elboi)
  }

# i <- 3; Yi <- data$Y[i, ]; Xi <- data$X[which(data$ij[, 1]==i), ]; mi <- eStep$M[i, ]; Si <- eStep$S[i, ]; xii <- eStep$xi[i, ]; beta <- true$beta; gamma <- true$gamma; C <- true$C ; logFactYi <- data$logFactY[i, ] 

ELBO <- function(data, mStep, eStep){
    elbo <- sum(sapply(1:nrow(data$Y), function(i){
      datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
      eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=eStep$S[i, ])
      ELBOi(datai=datai, mStep=mStep, eStepi=eStepi)
    }))
    return(elbo)
  }
  
ELBOold <- function(data, mStep, eStep){
    n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X); q <- ncol(eStep$M)
    nu <- matrix(data$X%*%mStep$gamma, n, p)
    mu <- matrix(data$X%*%mStep$beta, n, p)
    A <- exp(mu + eStep$M%*%t(mStep$C) +
               0.5*t(sapply(1:n, function(i){diag(mStep$C%*%diag(eStep$S[i, ])%*%t(mStep$C))})))
    elbo <- sum(nu*eStep$xi - log(1 + exp(nu)))
    elbo <- elbo - 0.5*(sum(eStep$M^2) + sum(eStep$S))
    elbo <- elbo + sum(eStep$xi * (-A + data$Y*(mu + eStep$M%*%t(mStep$C)) - data$logFactY))
    elbo <- elbo + sum(eStep$xi * log(eStep$xi/(1 - eStep$xi)) + log(1 - eStep$xi))
    elbo <- elbo + n*q/2 + 0.5*sum(log(eStep$S))
    return(elbo)
  }
  
################################################################################
# VE step
ElboMi <- function(mi, datai, mStep, eStepi){
  eStepi$mi <- mi
  ELBOi(datai=datai, mStep=mStep, eStepi=eStepi)
}
ElboSi <- function(Si, datai, mStep, eStepi){
  eStepi$Si <- Si
  ELBOi(datai=datai, mStep=mStep, eStepi=eStepi)
}
# ELBOMi <- function(mi, Yi, Xi, logFactYi, gamma, beta, C, xii, mi, Si){
#   ELBOi(mi=mi, Yi=Yi, X=Xi, logFactYi=logFactYi, gamma=gamma, beta=beta, C=C, xii=xii, Si=Si)
# }
# ElboMi <- function(mi, Yi, mui, Si, mStep, xii){
#   Ai <- exp(mui + mi%*%t(mStep$C) + 0.5*diag(mStep$C%*%diag(Si)%*%t(mStep$C)))
#   -0.5*sum(mi^2) -0.5*sum(Si) + sum(xii * (-Ai + Yi*(mui + mi%*%t(mStep$C))))
# }
# ElboSi <- function(Si, Yi, mui, mi, mStep, xii){ElboMi(mi=mi, Yi=Yi, mui=mui, Si=Si, mStep=mStep, xii=xii)}
ElboGradMi <- function(mi, datai, mStep, eStepi){
  mui <- as.vector(datai$Xi%*%mStep$beta)
  Ai <- exp(mui + mi%*%t(mStep$C) + 0.5*diag(mStep$C%*%diag(eStepi$Si)%*%t(mStep$C)))
  as.vector(-mi + (datai$Yi - Ai)%*%mStep$C)
}
ElboGradSi <- function(Si, datai, mStep, eStepi){
  mui <- as.vector(datai$Xi%*%mStep$beta)
  Ai <- exp(mui + eStepi$mi%*%t(mStep$C) + 0.5*diag(mStep$C%*%diag(Si)%*%t(mStep$C)))
  as.vector(0.5*1/Si - 1 - Ai%*%mStep$C)
}

VEstep <- function(data, mStep, eStep, tolXi=1e-4){
  n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X); q <- ncol(eStep$M)
  nu <- matrix(data$X%*%mStep$gamma, n, p)
  mu <- matrix(data$X%*%mStep$beta, n, p)
  A <- exp(mu + eStep$M%*%t(mStep$C) + 
             0.5*t(sapply(1:n, function(i){diag(mStep$C%*%diag(eStep$S[i, ])%*%t(mStep$C))})))
  xiLogit <- nu - A + data$Y*(eStep$M%*%t(mStep$C) + mu) - data$logFactY
  xi <- plogis(xiLogit)
  xi <- (xi + tolXi) / (1 + 2*tolXi)
  M <- t(sapply(1:n, function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=eStep$S[i, ])
    optim(par=eStep$M[i, ], fn=ElboMi, gr=ElboGradMi, 
          datai=datai, mStep=mStep, eStepi=eStepi, 
          method='BFGS', control=list(fnscale=-1))$par
    }))
  S <- t(sapply(1:n, function(i){
    datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], logFactYi=data$logFactY[i, ])
    eStepi <- list(xii=eStep$xi[i, ], mi=eStep$M[i, ], Si=eStep$S[i, ])
    optim(par=eStepi$Si, fn=ElboSi, gr=ElboGradSi, 
          datai=datai, mStep=mStep, eStepi=eStepi, 
          method='L-BFGS-B', control=list(fnscale=-1), lower=rep(1e-4, q))$par
    }))
  return(list(xi=xi, M=M, S=S))
}

################################################################################
# M step
ElboC <- function(vecC, data, mStep, eStep){
  mStep$C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboBeta <- function(beta, data, mStep, eStep){
  mStep$beta <- beta
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
ElboGamma <- function(gamma, data, mStep, eStep){
  mStep$gamma <- gamma
  ELBO(data=data, mStep=mStep, eStep=eStep)
}
# ElboC <- function(vecC, data, eStep, mu){
#   C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
#   A <- exp(mu + eStep$M%*%t(C) + 
#              0.5*t(sapply(1:n, function(i){diag(C%*%diag(eStep$S[i, ])%*%t(C))})))
#   sum(eStep$xi * (-A + data$Y*(mu + eStep$M%*%t(C))))
# }
ElboGradC <- function(vecC, data, mStep, eStep){
  C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
  mu <- matrix(data$X%*%mStep$beta, nrow(data$Y), ncol(data$Y))
  A <- exp(mu + eStep$M%*%t(C) + 
             0.5*t(sapply(1:n, function(i){diag(C%*%diag(eStep$S[i, ])%*%t(C))})))
  as.vector(t(eStep$xi*(data$Y - A))%*%eStep$M - (t(A)%*%eStep$S)*C)
}
# ElboBeta <- function(beta, data, eStep, C){
#   mu <- matrix(data$X%*%beta, nrow(data$Y), ncol(data$Y))
#   A <- exp(mu + eStep$M%*%t(C) + 
#              0.5*t(sapply(1:n, function(i){diag(C%*%diag(eStep$S[i, ])%*%t(C))})))
#   sum(eStep$xi * (data$Y*(mu + eStep$M%*%t(C)) - A))
# }
ElboGradBeta <- function(beta, data, mStep, eStep){
  mu <- matrix(data$X%*%beta, nrow(data$Y), ncol(data$Y))
  A <- exp(mu + eStep$M%*%t(mStep$C) + 
             0.5*t(sapply(1:n, function(i){diag(mStep$C%*%diag(eStep$S[i, ])%*%t(mStep$C))})))
  as.vector(t(data$X) %*% as.vector(data$Y - A))
}
# ElboGamma <- function(gamma, data, eStep){
#   nu <- matrix(data$X%*%gamma, nrow(data$Y), ncol(data$Y))
#   sum(eStep$xi * nu - log(1 + exp(nu)))
# }
ElboGradGamma <- function(gamma, data, mStep, eStep){
  nu <- matrix(data$X%*%gamma, nrow(data$Y), ncol(data$Y))
  probU <- plogis(nu)
  as.vector(t(data$X)%*%as.vector(eStep$xi - probU))
}

Mstep <- function(data, mStep, eStep){
  mu <- matrix(data$X%*%mStep$beta, n, p)
  vecC <- optim(par=as.vector(mStep$C), fn=ElboC, gr=ElboGradC, data=data, mStep=mStep, eStep=eStep, 
        method='BFGS', control=list(fnscale=-1))$par
  C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
  beta <- optim(par=mStep$beta, fn=ElboBeta, gr=ElboGradBeta, data=data, mStep=mStep, eStep=eStep, 
                method='BFGS', control=list(fnscale=-1))$par
  gamma <- optim(par=mStep$gamma, fn=ElboGamma, gr=ElboGradGamma, data=data, mStep=mStep, eStep=eStep, 
                method='BFGS', control=list(fnscale=-1))$par
  return(list(gamma=gamma, beta=beta, C=C))
}
