SimZiPLNmiss <- function(sim, obs=1){
  sim$Omega <- matrix(rbinom(prod(dim(sim$Y)), 1, obs), nrow(sim$Y), ncol(sim$Y))
  sim$Yfull <- sim$Y
  sim$Y <- sim$Y * sim$Omega
  return(sim)
}

InitZiPLNold <- function(data, q){
  pres <- which(data$Omega==1)
  reg <- lm(as.vector(log(1+data$Y)[pres]) ~ -1 + data$X[pres, ])
  res <- matrix(0, nrow(data$Y), ncol(data$Y))
  res[pres] <- reg$residuals
  pca <- prcomp(res, rank=q)
  mStep <- list(gamma=as.vector(glm(as.vector(1*(data$Y > 0)) ~ -1 + data$X, family='binomial')$coef), 
                beta=reg$coef, 
                C=pca$rotation %*% diag(pca$sdev[1:q, drop=FALSE]))
  eStep <- list(xi=matrix(sum(data$Omega*(data$Y > 0))/sum(data$Omega), n, ncol(data$Y)), M=matrix(0, n, q), S=matrix(1e-4, n, q))
  return(list(mStep=mStep, eStep=eStep, reg=reg, pca=pca))
}
