# Tests for the ZIP distrib

library('bizicount')
source('../Functions/FunctionsZIP.R')

# Parms 
pi <- 0.1; lambda <- 5

# Sim 
B <- 1e4
Y <- rzip(n=B, psi=1-pi, lambda=lambda)
c(mean(Y), pi*lambda)
c(var(Y), pi*lambda + pi*(1-pi)*lambda^2)
qzip(p=c(.05, .95), psi=1-pi, lambda=lambda)
quantile(p=c(.05, .95), Y)
