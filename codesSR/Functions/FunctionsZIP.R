# EM algo for ZIP distribution with covariates
#########################################################################
# !!! the logistic part codes for ABSENCE => gamma <- -gamma in the end
#########################################################################

# Initialisation
initZIP <- function(X, Y, Y0){
  gamma <- glm(Y0 ~ -1 + X, family='binomial')$coef
  beta <- glm(Y ~ -1 + X, family='poisson')$coef
  return(list(gamma=gamma, beta=beta))
}
# Fonctions
objGamma <- function(gamma, Y0, X, tau){
  (t(tau)%*%X%*%gamma - sum(log(1+exp(X%*%gamma))))[1, 1]
}
gradGamma <- function(gamma, Y0, X, tau){
  as.vector(t(X)%*%(tau-plogis(X%*%gamma)))
}
objBeta <- function(beta, Y, X, tau){
  (-sum((1-tau)*exp(X%*%beta)) + ((1-tau)*Y)%*%X%*%beta)[1, 1]
}
gradBeta <- function(beta, Y, X, tau){
  as.vector(t(X)%*%((1-tau)*(Y - exp(X%*%beta))))
}
# logLik <- function(gamma, beta, Y, Y0, X){
#   pi <- plogis(X%*%gamma); lambda <- exp(X%*%beta)
#   sum(log(pi*Y0 + (1-pi)*dpois(Y, lambda)))
# }
logLik <- function(gamma, beta, Y, Y0, X){
  pres <- which(Y0==0)
  pi <- plogis(X%*%gamma); lambda <- exp(X%*%beta)
  logL <- log(pi*Y0 + (1-pi)*dpois(Y, lambda))
  logL[pres] <- log(1-pi[pres])+dpois(Y[pres], lambda[pres], log=TRUE)
  sum(logL)
}

# Algorithme EM
EmZIP <- function(X, Y, tol=1e-4, iterMax=1e3, plot=FALSE){
  Y0 <- 1*(Y==0)
  start <- initZIP(X, Y, Y0)
  gamma <- start$gamma; beta <- start$beta
  diff <- 2*tol; iter <- 0
  logL <- rep(NA, iterMax)
  while((diff > tol) & (iter < iterMax)){
    iter <- iter+1
    # Etape E
    pi <- plogis(X%*%gamma); lambda <- exp(X%*%beta)
    tauNew <- as.vector(pi*Y0 / (pi + (1-pi)*exp(-lambda)))
    # Etape M
    gammaNew <- optim(par=gamma, fn=objGamma, gr=gradGamma, Y0=Y0, X=X, tau=tauNew, 
                      control=list(fnscale=-1))$par
    betaNew <- optim(par=beta, fn=objBeta, gr=gradBeta, Y=Y, X=X, tau=tauNew, 
                     control=list(fnscale=-1))$par
    # Test & mise à jour
    if(iter > 1){diff <- max(max(abs(tauNew-tau)), max(gammaNew-gamma), max(betaNew-beta))}
    logL[iter] <- logLik(gamma, beta, Y, Y0, X)
    gamma <- gammaNew; beta <- betaNew; tau <- tauNew
    if(plot){if(iter%%round(sqrt(iterMax))==0){
      plot(logL[1:iter], ylim=quantile(logL[1:iter], probs=c(0.1, 1)), type='b', main='ZIP')
    }}
  }
  if(plot){
    plot(logL[1:iter], ylim=quantile(logL[1:iter], probs=c(0.1, 1)), type='b', main='ZIP')
  }
  return(list(gamma=-gamma, beta=beta, iter=iter, logLpath=logL[1:iter], logL=logL[iter]))
}
