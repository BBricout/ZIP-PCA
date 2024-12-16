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
Psi2Estep <- function(psi, n, d, p, q){
  list(M=matrix(psi[(1:(n*q))], n, q, byrow=TRUE), 
       S=matrix(psi[(n*q) + (1:(n*q))], n, q, byrow=TRUE))
}
ComputeXi <- function(data, mStep, eStep, tolXi=1e-4){
  nuMuA <- NuMuA(data=data, mStep=mStep, eStep=eStep)
  xi <- plogis(nuMuA$nu - data$Omega*nuMuA$A)
  xi <- (xi + tolXi) / (1 + 2*tolXi)
  xi[which(data$Omega*data$Y>0)] <- 1
  return(xi)
}
MakeCOrtho <- function(C){
  q <- ncol(C); eigC <- eigen(C%*%t(C))
  return(eigC$vectors[, 1:q]%*%diag(sqrt(eigC$values[1:q])))
}

################################################################################
# Simul
SimZiPLN <- function(n, p, d, q, coefC=1, beta=NULL, gamma=NULL, C=NULL, X=NULL){
  if(is.null(X)){X <- matrix(rnorm(n*p*d), n*p, d); X[, 1] <- 1}
  ij <- cbind(rep(1:n, p), rep(1:p, each=n))
  # presence
  if(is.null(gamma)){gamma <- rnorm(d)/sqrt(d)}
  nu <- matrix(X%*%gamma, n, p)
  probU <- plogis(nu)
  U <- matrix(rbinom(n*p, 1, probU), n, p)
  # latent
  W <- matrix(rnorm(n*q), n, q)
  if(is.null(C)){C <- matrix(rnorm(p*q), p, q)/sqrt(q)}
  Z <- W%*%t(C)
  # abundance
  if(is.null(beta)){beta <- rnorm(d)/sqrt(d)}
  mu <- matrix(X%*%beta, n, p)
  lambdaY <- exp(mu + Z)
  Yall <- matrix(rpois(n*p, lambdaY), n, p)
  Y <- Yall*U
  return(list(X=X, Y=Y, ij=ij, U=U, W=W, Z=Z, Yall=Yall, gamma=gamma, beta=beta, C=C))
}

