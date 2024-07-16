# seed <- .Random.seed
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIPLNmiss.R')

Parms2Steps <- function(parms, data){
  n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X)
  mStep <- list(gamma=parms[1:d], beta=parms[d+(1:d)], 
                C=matrix(parms[(2*d)+(1:(p*q))], p, q))
  eStep <- list(M=matrix(parms[(2*d)+(p*q)+(1:(n*q))], n, q, byrow=TRUE), 
                S=matrix(parms[(2*d)+(p*q)+(n*q)+(1:(n*q))], n, q, byrow=TRUE))
  eStep$xi <- ComputeXi(data=data, mStep=mStep, eStep=eStep, tolXi=0)
  return(list(m=mStep, e=eStep))  
}


Parms2Params <- function(parms, data){
  n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X)
  steps <- Parms2Steps(parms=parms, data=data)
  return(list(B=as.matrix(steps$m$beta), 
              D=as.matrix(steps$m$gamma), 
              C=as.matrix(steps$m$C), 
              M=as.matrix(steps$e$M), 
              S=as.matrix(steps$e$S)))
}
















