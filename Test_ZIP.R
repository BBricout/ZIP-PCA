library(Rcpp)
library(Matrix)
library(mvtnorm)
library(stringr)
library(ggplot2)
library(matrixcalc)
library(missForest)

local.dir = getwd()
sourceCpp(file.path(local.dir,"src/optim_rank_ZIP.cpp"))

#------------ Simulation d'un jeu de données -----------------


n <- 500
p <- 30
d <- 2
q <- 3

X <- cbind(c(rep(1, n*p)), matrix(rnorm(n*p*d), nrow = n*p))
B <- c(2, 1, 0.5)
D <- c(-0.5, 0.5, -0.7)

mu <- VectorToMatrix(X%*%B, n, p)
nu <- VectorToMatrix(X%*%D, n, p)
 
Prob <- plogis(nu)
U <- matrix(rbinom(n*p, p = Prob, size = 1), nrow = n)

W <- matrix(rnorm(n*q), nrow = n)/2
C <- matrix(rnorm(p*q,0,1), nrow = p)/ sqrt(q)

var(X%*%B)
var(MatrixToVector(W%*%t(C)))

Lambda <- exp(mu + W %*% t(C))
Z <- matrix(rpois(n*p, lambda = Lambda), nrow = n)

Y <- ifelse(U == 0, 0, Z)
range(Y)



#-------------- Test de l'optim ----------------

Init_ZIP <- function(Y, X, q){
  
  n <- nrow(Y)
  p <- ncol(Y)
  vecY <- MatrixToVector(Y)
  
  U <- ifelse(Y == 0, 0, 1)
  fit.logit <- glm(vec(U) ~ -1 + X, family = "binomial", na.action = na.exclude)
  D <- as.matrix(fit.logit$coefficients)
  
  fit <- lm(log(1 + vecY) ~ -1 + X, na.action = na.exclude)
  B <- as.matrix(fit$coefficients)
  res.mat <- VectorToMatrix(fit$residuals, n, p)
  
  svdM <- svd(res.mat, nu = q, nv = p)
  
  C <- svdM$v[, 1:q, drop = FALSE] %*% diag(svdM$d[1:q], nrow = q, ncol = q)/sqrt(n)
  M  <- svdM$u[, 1:q, drop = FALSE] %*% diag(svdM$d[1:q], nrow = q, ncol = q) %*% t(svdM$v[1:q, 1:q, drop = FALSE])
  S <- matrix(0.1, n, q)
  
  return(list(B = B, D = D, C = C, M = M, S = S))
}



data <- list(Y = Y, R = matrix(1, n, p), X = X)
params <- Init_ZIP(Y, X, q)
config <- PLNPCA_param()$config_optim


out <- nlopt_optimize_ZIP(data, params, config)

elbo <- out$objective_values
plot(elbo, type='b', xlab='iter', ylim=quantile(elbo, probs=c(0.1, 0.9)))
length(out$objective_values)

sigma <- C%*%t(C)
sigma.out <- out$C%*%t(out$C)

plot(X%*%B, X%*%out$B, main = "B") ; abline(0,1)
plot(X%*%D, X%*%out$D, main = "D") ; abline(0,1)
plot(sigma, sigma.out, main = "sigma") ; abline(0,1)


M <- params$M
S <- params$S

A <- exp(mu + M%*%t(C) + 0.5 * S %*% t(C * C))

E <- ifelse(Y==0, nu - A, 1)
range(mu)
mu.init <- VectorToMatrix(X%*%params$B, n, p)
range(mu.init)
range(exp(mu))
range(exp(-mu))












