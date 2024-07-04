library(Rcpp)
library(Matrix)
library(mvtnorm)
library(stringr)
library(ggplot2)
library(matrixcalc)
library(missForest)
library(PLNmodels)
library(patchwork)

local.dir = getwd()
sourceCpp(file.path(local.dir,"src/optim_rank_ZIP.cpp"))
sourceCpp("src/optim_rank_cov_inter.cpp")
source("Initialisations.R")

Y <- as.matrix(read.csv("data/Colvert_count.csv"))
Y <- Y[,-1]
head(Y)

X <- as.matrix(read.csv("data/Colvert_covariate.csv"))
X <- X[,-1]
head(X)


range(Y)

n <- nrow(Y)
p <- ncol(Y)
q <- 6

# Initialisations

config <- PLNPCA_param()$config_optim


params <- lapply(1:p, function(i)
  Init_PLNPCA(Y, X, i))

params.ZIP <- lapply(1:p, function(i)
  Init_ZIP(Y, X, i))


O <- matrix(0, nrow = n, ncol = p)  # n*p
w <- rep(1, n)

data <- list(Y = Y, R = matrix(1, n, p), X = X, O = O, w = w)

data.ZIP <- list(Y = Y, R = matrix(1, n, p), X = X)

# out

out <- lapply(1:p, function(i)
  nlopt_optimize_rank_cov(data, params[[i]], config))

out.ZIP <- lapply(1:p, function(i)
  nlopt_optimize_ZIP(data.ZIP, params.ZIP[[i]], config))

# Elbo

Elbo <- lapply(1:p, function(j)
  sum(out[[j]]$Ji))

plot(unlist(Elbo))

Elbo.ZIP <- lapply(1:p, function(j)
  sum(out[[j]]$objective_values[[length(out[[j]]$objective_values)]]))

plot(unlist(Elbo.ZIP))

# Predictions

XB <- lapply(1:p, function(i)
  VectorToMatrix(X%*%out[[i]]$B, n, p))
A.pred <- lapply(1:p, function(i)
  XB[[i]] + out[[i]]$M %*% t(out[[i]]$C) + 1/2 * hadamard.prod(out[[i]]$S, out[[i]]$S) %*% t(hadamard.prod(out[[i]]$C, out[[i]]$C)))
pred <- lapply(1:p, function(i)
  exp(A.pred[[i]]))

pred.ZIP <- lapply(1:p, function(i)
  out.ZIP[[i]]$xi * out.ZIP[[i]]$A)


data.frame(
  fitted   = as.vector(1 + pred[[7]]),
  observed = as.vector(1 + Y)
) %>%
  ggplot(aes(x = observed, y = fitted)) +
  geom_point(size = .5, alpha =.25 ) +
  scale_x_log10(limits = c(1,NA)) +
  scale_y_log10(limits = c(1,NA)) +
  geom_abline(slope = 1, intercept = 0, col ="grey70", linetype = "dashed") +
  theme_bw() + annotation_logticks() +
  ggtitle("Prédictions modèle PLNPCA")-> plot.ref

data.frame(
  fitted   = as.vector(1 + pred.ZIP[[7]]),
  observed = as.vector(1 + Y)
) %>%
  ggplot(aes(x = observed, y = fitted)) +
  geom_point(size = .5, alpha =.25 ) +
  scale_x_log10(limits = c(1,NA)) +
  scale_y_log10(limits = c(1,NA)) +
  geom_abline(slope = 1, intercept = 0, col ="grey70", linetype = "dashed") +
  theme_bw() + annotation_logticks() +
  ggtitle("Prédictions modèle ZIP PLNPCA")-> plot.ZIP


(plot.ref + plot.ZIP & theme(legend.position = "bottom")) + plot_layout(guides = "collect")








