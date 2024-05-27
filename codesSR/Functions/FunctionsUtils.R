SimZiPLN <- function(n, p, d, q, beta0=2){
  X <- matrix(rnorm(n*p*d), n*p, d); X[, 1] <- 1
  ij <- cbind(rep(1:n, p), rep(1:p, each=n))
  # presence
  gamma <- rnorm(d)/sqrt(d)
  nu <- matrix(X%*%gamma, n, p)
  probU <- plogis(nu)
  U <- matrix(rbinom(n*p, 1, probU), n, p)
  # latent
  W <- matrix(rnorm(n*q), n, q)
  C <- matrix(rnorm(p*q), p, q)/sqrt(q)
  Z <- W%*%t(C)
  # abundance
  beta <- rnorm(d)/sqrt(d); beta[1] <- beta[1] + beta0
  mu <- matrix(X%*%beta, n, p)
  lambdaY <- exp(mu + Z)
  Yall <- matrix(rpois(n*p, lambdaY), n, p)
  Y <- Yall*U
  return(list(X=X, Y=Y, ij=ij, U=U, W=W, Z=Z, Yall=Yall, gamma=gamma, beta=beta, C=C))
}

