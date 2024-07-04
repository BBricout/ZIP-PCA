#-------------Fonctions-------------

library(Rcpp)
library(Matrix)
library(mvtnorm)
library(stringr)
library(ggplot2)
library(matrixcalc)
library(missMDA)
library(PLNmodels)

sourceCpp("src/optim_rank_ZIP.cpp")

#---------Initialisations-------------

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


#---------------Fonctions-------------


Miss.ZIPPCA <- function(Y, # Table de comptages n*p qui peut contenir des données manquantes
                        X, # Covariables np*d dont une colonne de 1 pour l'intercept
                        q, # Dimension de l'espace latent q
                        params = NULL, # Paramètres fourni en entrée
                        config = NULL,
                        tolXi = NULL){
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  if (is.null(params)){params <- Init_ZIP(Y, X, q)}
  if (is.null(config)){config <- PLNPCA_param()$config_optim}
  if (is.null(tolXi)){tolXi <- 1e-4}
  
  R <- ifelse(is.na(Y), 0, 1) # Masque qui met des 0 à la place des données manquantes
  
  Y.na <- ifelse(R == 0, 0, Y)
  
  data <- list(Y = Y.na,
               R = R,
               X = X)
  
  out <- nlopt_optimize_ZIP(data, params, config)
  
  mu <- VectorToMatrix(X%*%out$B, n, p)
  nu <- VectorToMatrix(X%*%out$D, n, p)
  
  mStep <- list(gamma = out$D, beta = out$B, C = out$C)
  eStep <- list(M = out$M, S = out$S,  xi = out$xi)
  pred <- list(A = out$A, nu = nu, mu = mu)
  iter <- out$monitoring$iterations
  elboPath <- out$objective_values
  elbo <- out$objective_values[length(out$objective_values)]

  
  res <- list(mStep = mStep, 
              eStep = eStep, 
              pred = pred, 
              iter = iter, 
              elboPath = elboPath, 
              elbo = elbo,
              params.init = params,
              monitoring = out$monitoring)
  
  return(res)
  
}










