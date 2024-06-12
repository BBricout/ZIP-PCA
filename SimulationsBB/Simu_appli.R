library(Rcpp)
library(Matrix)
library(mvtnorm)
library(stringr)
library(ggplot2)
library(matrixcalc)
library(missForest)
library(PLNmodels)

local.dir = getwd()
sourceCpp(file.path(local.dir,"src/optim_rank_ZIP.cpp"))

load("Simulations/Simu_data.Rdata")

Y <- Data$Y
Z <- Data$Z
C <- Data$C
B <- Data$B
X <- Data$X
W <- Data$W
mu <- Data$mu
nu <- Data$nu
D <- Data$D
Y.na10 <- Data$Y.na10
Y.na30 <- Data$Y.na30
Y.na60 <- Data$Y.na60

n <- 300
p <- 30
d <- 2
q <- 3

config <- PLNPCA_param()$config_optim

#------------- Initialisation -----------------




# Paramètres initiaux

params <- Init_ZIP(Y, X, q)
params.na10 <- Init_ZIP(Y.na10, X, q)
params.na30 <- Init_ZIP(Y.na30, X, q)
params.na60 <- Init_ZIP(Y.na60, X, q)

# Données

data <- list(Y = Y, R = matrix(1, n, p), X = X)
data.na10 <- list(Y = Y, R = ifelse(is.na(Y.na10), 0, 1), X = X)
data.na30 <- list(Y = Y, R = ifelse(is.na(Y.na30), 0, 1), X = X)
data.na60 <- list(Y = Y, R = ifelse(is.na(Y.na60), 0, 1), X = X)


# Ajustement

out <- nlopt_optimize_ZIP(data, params, config)
out.na10 <- nlopt_optimize_ZIP(data.na10, params.na10, config)
out.na30 <- nlopt_optimize_ZIP(data.na30, params.na30, config)
out.na60 <- nlopt_optimize_ZIP(data.na60, params.na60, config)

# Visualisation des regresseurs

par(mfrow = c(2, 2))  # 2 lignes, 2 colonnes

plot(B, out$B, main = "B en données complètes") ; abline(0,1)
plot(B, out.na10$B, main = "B pour 10% de NA") ; abline(0,1)
plot(B, out.na30$B, main = "B pour 30% de NA") ; abline(0,1)
plot(B, out.na60$B, main = "B pour 60% de NA") ; abline(0,1)


# W%*%t(C)

Z.ref <- W%*%t(C)
Z.out <- out$M%*%t(out$C)
Z.out.na10 <- out.na10$M%*%t(out.na10$C)
Z.out.na30 <- out.na30$M%*%t(out.na30$C)
Z.out.na60 <- out.na60$M%*%t(out.na60$C)

par(mfrow = c(2, 2))
plot(Z.ref, Z.out, main = "Z") ; abline(0,1)
plot(Z.ref, Z.out.na10, main = "Z 10% de NA") ; abline(0,1)
plot(Z.ref, Z.out.na30, main = "Z 30% de NA") ; abline(0,1)
plot(Z.ref, Z.out.na60, main = "Z 60% de NA") ; abline(0,1)

# Predictions

par(mfrow = c(2, 2))
plot(log(1 + Y), log(1 + out$xi * out$A), main = "Prédictions complètes") ; abline(0,1)
plot(log(1 + Y), log(1 + out.na10$xi * out.na10$A), main = "Prédictions 10% de NA") ; abline(0,1)
plot(log(1 + Y), log(1 + out.na30$xi * out.na30$A), main = "Prédictions 30% de NA") ; abline(0,1)
plot(log(1 + Y), log(1 + out.na60$xi * out.na60$A), main = "Prédictions 60% de NA") ; abline(0,1)












