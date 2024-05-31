# Tests for the ZIP distrib


################################################################################
# Regression
library(pscl)
source('../Functions/FunctionsZIP.R')

n <- 1e3; d <- 3; 
X <- matrix(rnorm(n*d), n, d); X[, 1] <- 1
gamma <- rnorm(d); beta <- rnorm(d); beta[1] <- beta[1] + 2
lambda <- exp(X%*%beta); Y <- rpois(n, lambda)
pi <- plogis(X%*%gamma); Z <- rbinom(n, 1, pi)
Y <- Y*(1-Z)
zip1 <- zeroinfl(Y ~ -1 + X, dist='poisson')
zip2 <- EmZIP(X=X, Y=Y)
# rbind(gamma, -zip1$coefficients$zero, -zip2$gamma)
# rbind(beta, zip1$coefficients$count, zip2$beta)
plot(-zip1$coefficients$zero, zip2$gamma); abline(0, 1)
plot(zip1$coefficients$count, zip2$beta); abline(0, 1)
c(zip1$loglik, zip2$logL)

# ################################################################################
# # Quantile
# library(bizicount); 
# # Parms 
# pi <- 0.1; lambda <- 5
# 
# # Sim 
# B <- 1e4
# Y <- rzip(n=B, psi=1-pi, lambda=lambda)
# c(mean(Y), pi*lambda)
# c(var(Y), pi*lambda + pi*(1-pi)*lambda^2)
# qzip(p=c(.05, .95), psi=1-pi, lambda=lambda)
# quantile(p=c(.05, .95), Y)
