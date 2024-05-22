# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
# seed <- 1; set.seed(seed)
source('FunctionsZI-PLN.R')

# Parms
n <- 500; d <- 10; p <-30; q <- 3
sim <- SimZiPLN(n=n, p=p, d=d, q=q)
data <- list(X=sim$X, Y=sim$Y, ij=sim$ij, logFactY=lgamma(sim$Y+1))
true <- list(gamma=sim$gamma, beta=sim$beta, C=sim$C, 
             U=sim$U, W=sim$W, Z=sim$Z,
             M=sim$W, S=matrix(1e-4, n, q), xi=matrix(plogis(sim$X%*%sim$gamma), n, p))

# Init
# mStep <- list(gamma=true$gamma, beta=true$beta, C=true$C)
# eStep <- list(xi=true$xi, M=true$M, S=true$S)
init <- InitZiPLN(data)
mStep <- init$mStep; eStep <- init$eStep
ELBO(data=data, eStep=eStep, mStep=mStep)
ELBOold(data=data, eStep=eStep, mStep=mStep)

# Fit
iterMax <- 1e3; tol <- 1e-6
elboPath <- rep(NA, iterMax)
diff <- 2*tol; iter <- 1
elboPath[iter] <- ELBO(data=data, mStep=mStep, eStep=eStep)
while((iter < iterMax) & (diff > tol)){
  cat(iter, elboPath[iter], diff, '\n')
  iter <- iter+1
  Mold <- eStep$M; Cold <- mStep$C
  mStep <- Mstep(data=data, mStep=mStep, eStep=eStep)
  # mStep$gamma <- true$gamma; mStep$beta <- true$beta; mStep$C <- true$C 
  # mStep <- list(gamma=true$gamma, beta=true$beta, C=true$C)
  # eStep <- VEstep(data=data, mStep=mStep, eStep=eStep)
  # eStep <- list(xi=true$xi, M=true$M, S=true$S)
  eStep$xi <- true$xi; eStep$M <- true$M; eStep$S <- true$S
  elboPath[iter] <- ELBO(data=data, mStep=mStep, eStep=eStep)
  # diff <- max(abs(eStep$M - Mold)); Mold <- eStep$M
  diff <- max(abs(mStep$C - Cold)); Cold <- mStep$C
  if(iter%%round(sqrt(iterMax))==0){
    plot(elboPath[1:iter], type='b', ylim=quantile(elboPath[1:iter], probs=c(0.1, 1), na.rm=TRUE))
  }
}
cat(iter, elboPath[iter], diff, '\n')
elboPath <- elboPath[1:iter]

par(mfrow=c(3, 2))
plot(elboPath, type='b', ylim=quantile(elboPath, probs=c(0.1, 1)))
plot(true$gamma, mStep$gamma); abline(a=0, b=1, h=0, v=0)
plot(true$beta, mStep$beta); abline(a=0, b=1, h=0, v=0)
plot(true$C, mStep$C); abline(a=0, b=1, h=0, v=0)
plot(true$C%*%t(true$C), mStep$C%*%t(mStep$C)); abline(a=0, b=1, h=0, v=0)
plot(true$Z, eStep$M%*%t(mStep$C)); abline(a=0, b=1, h=0, v=0)
