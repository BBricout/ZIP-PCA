Bobj <- function(parms, data, tolXi){
  params <- Parms2Params(parms=parms, data=data, tolXi)
  Elbo(data=data, params=params, tolXi)$objective
}
Bobj.neg <- function(parms, data, tolXi){
  -Bobj(parms, data, tolXi)
}
Bgrad <- function(parms, data, tolXi){
  params <- Parms2Params(parms=parms, data=data, tolXi)
  grads <- Grad(data=data, params=params, tolXi)
  c(as.vector(grads$gradD), as.vector(grads$gradB), as.vector(grads$gradC), 
    as.vector(t(grads$gradM)), as.vector(t(grads$gradS)))
}
Bgrad.neg <- function(parms, data, tolXi){
  -Bgrad(parms, data, tolXi)
}
BobjLogS <- function(parmsLogS, data, tolXi){
  params <- Parms2Params(parms=parmsLogS, data=data,  tolXi)
  params$logS <- params$S
  Elbo_grad_logS_Rcpp(data, params, tolXi)$objective
}
BobjLogS.neg <- function(parmsLogS, data, tolXi){
  -BobjLogS(parmsLogS, data, tolXi)
}
BgradLogS <- function(parmsLogS, data, tolXi){
  params <- Parms2Params(parms=parmsLogS, data=data, tolXi)
  params$logS <- params$S
  elbo <- Elbo_grad_logS_Rcpp(data, params, tolXi)
  c(as.vector(elbo$gradD), as.vector(elbo$gradB), as.vector(elbo$gradC),
    as.vector(t(elbo$gradM)), as.vector(t(elbo$gradS))*as.vector(t(params$S)))
}
BgradLogS.neg <- function(parmsLogS, data, tolXi){
  -BgradLogS(parmsLogS, data, tolXi)
}

