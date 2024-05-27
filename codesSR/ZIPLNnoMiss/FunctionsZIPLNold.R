# Functions ZI-PLN

ELBOold <- function(data, mStep, eStep){
    n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X); q <- ncol(eStep$M)
    nu <- matrix(data$X%*%mStep$gamma, n, p)
    mu <- matrix(data$X%*%mStep$beta, n, p)
    A <- exp(mu + eStep$M%*%t(mStep$C) +
               0.5*t(sapply(1:n, function(i){diag(mStep$C%*%diag(eStep$S[i, ])%*%t(mStep$C))})))
    elbo <- sum(nu*eStep$xi - log(1 + exp(nu)))
    elbo <- elbo - 0.5*(sum(eStep$M^2) + sum(eStep$S))
    elbo <- elbo + sum(eStep$xi * (-A + data$Y*(mu + eStep$M%*%t(mStep$C)) - data$logFactY))
    elbo <- elbo + sum(eStep$xi * log(eStep$xi/(1 - eStep$xi)) + log(1 - eStep$xi))
    elbo <- elbo + n*q/2 + 0.5*sum(log(eStep$S))
    return(elbo)
  }
  
################################################################################
# VE step
# ELBOMi <- function(mi, Yi, Xi, logFactYi, gamma, beta, C, xii, mi, Si){
#   ELBOi(mi=mi, Yi=Yi, X=Xi, logFactYi=logFactYi, gamma=gamma, beta=beta, C=C, xii=xii, Si=Si)
# }
# ElboMi <- function(mi, Yi, mui, Si, mStep, xii){
#   Ai <- exp(mui + mi%*%t(mStep$C) + 0.5*diag(mStep$C%*%diag(Si)%*%t(mStep$C)))
#   -0.5*sum(mi^2) -0.5*sum(Si) + sum(xii * (-Ai + Yi*(mui + mi%*%t(mStep$C))))
# }
# ElboSi <- function(Si, Yi, mui, mi, mStep, xii){ElboMi(mi=mi, Yi=Yi, mui=mui, Si=Si, mStep=mStep, xii=xii)}

################################################################################
# M step
# ElboC <- function(vecC, data, eStep, mu){
#   C <- matrix(vecC, ncol(data$Y), ncol(eStep$M))
#   A <- exp(mu + eStep$M%*%t(C) + 
#              0.5*t(sapply(1:n, function(i){diag(C%*%diag(eStep$S[i, ])%*%t(C))})))
#   sum(eStep$xi * (-A + data$Y*(mu + eStep$M%*%t(C))))
# }
# ElboBeta <- function(beta, data, eStep, C){
#   mu <- matrix(data$X%*%beta, nrow(data$Y), ncol(data$Y))
#   A <- exp(mu + eStep$M%*%t(C) + 
#              0.5*t(sapply(1:n, function(i){diag(C%*%diag(eStep$S[i, ])%*%t(C))})))
#   sum(eStep$xi * (data$Y*(mu + eStep$M%*%t(C)) - A))
# }
# ElboGamma <- function(gamma, data, eStep){
#   nu <- matrix(data$X%*%gamma, nrow(data$Y), ncol(data$Y))
#   sum(eStep$xi * nu - log(1 + exp(nu)))
# }
