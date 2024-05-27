SimZiPLNmiss <- function(sim, obs=1){
  sim$Omega <- matrix(rbinom(prod(dim(sim$Y)), 1, obs), nrow(sim$Y), ncol(sim$Y))
  sim$Yfull <- sim$Y
  sim$Y <- sim$Y * sim$Omega
  return(sim)
}
