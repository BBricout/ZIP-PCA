# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
# seed <- 1; set.seed(seed)
source('FunctionsZI-PLN.R')
library(PLNmodels)

# Parms
n <- 100; d <- 5; p <- 10; q <- 2
sim <- SimZiPLN(n=n, p=p, d=d, q=q)
data <- list(X=sim$X, Y=sim$Y, ij=sim$ij, logFactY=lgamma(sim$Y+1))
true <- list(gamma=sim$gamma, beta=sim$beta, C=sim$C, 
             U=sim$U, W=sim$W, Z=sim$Z,
             M=sim$W, S=matrix(1e-4, n, q), xi=matrix(plogis(sim$X%*%sim$gamma), n, p))
c(sum(diag(cov(matrix(sim$X%*%sim$beta, n, p)))), sum(diag(cov(sim$Z))))

# Oracle
oracle <- OracleZiPLN(sim)

# Init
mStep <- eStep <- true
init <- InitZiPLN(data); mStep <- init$mStep; eStep <- init$eStep
# pln <- PLN(as.vector(data$Y) ~ -1 + data$X); eStep$beta <- as.vector(pln$model_par$B)

# # First iterations
# ELBO(data=data, mStep=mStep, eStep=eStep)
# ElboGradC(vecC=as.vector(mStep$C), data=data, mStep=mStep, eStep=eStep)
# mStep <- Mstep(data=data, mStep=mStep, eStep=eStep)
# ELBO(data=data, mStep=mStep, eStep=eStep)
# eStep <- VEstep(data=data, mStep=mStep, eStep=eStep)
# ELBO(data=data, mStep=mStep, eStep=eStep)

# Fit
iterMax <- 1e3; tol <- 1e-4
elboPath <- rep(NA, iterMax)
diff <- 2*tol; iter <- 1
elboPath[iter] <- ELBO(data=data, mStep=mStep, eStep=eStep)
cat(iter, elboPath[iter], diff, '\n')
while((iter < iterMax) & (diff > tol)){
  iter <- iter+1
  Mold <- eStep$M; Cold <- mStep$C
  mStep <- Mstep(data=data, mStep=mStep, eStep=eStep)
  # mStep$gamma <- true$gamma; mStep$beta <- true$beta; # mStep$C <- true$C
  # mStep <- list(gamma=true$gamma, beta=true$beta, C=true$C)
  eStep <- VEstep(data=data, mStep=mStep, eStep=eStep)
  # eStep <- list(xi=true$xi, M=true$M, S=true$S)
  # eStep$xi <- true$xi; eStep$M <- true$M; eStep$S <- true$S
  elboPath[iter] <- ELBO(data=data, mStep=mStep, eStep=eStep)
  cat(iter, 'E', elboPath[iter], diff, '\n')
  diffM <- max(abs(eStep$M - Mold)); Mold <- eStep$M
  diffC <- max(abs(mStep$C - Cold)); Cold <- mStep$C
  diff <- max(diffM, diffC)
  if(iter%%round(sqrt(iterMax))==0){
    plot(elboPath[1:iter], type='b', ylim=quantile(elboPath[1:iter], probs=c(0.1, 1), na.rm=TRUE))
  }
}
elboPath <- elboPath[1:iter]

par(mfrow=c(3, 3))
plot(elboPath, type='b', ylim=quantile(elboPath, probs=c(0.1, 1)))
plot(true$gamma, mStep$gamma); abline(a=0, b=1, h=0, v=0)
points(true$gamma, oracle$mStep$gamma, col=2)
plot(true$beta, mStep$beta); abline(a=0, b=1, h=0, v=0)
points(true$beta, oracle$mStep$beta, col=2)
plot(true$C%*%t(true$C), mStep$C%*%t(mStep$C)); abline(a=0, b=1, h=0, v=0)
plot(true$C%*%t(true$C), oracle$mStep$C%*%t(oracle$mStep$C), col=2); abline(a=0, b=1, h=0, v=0)
plot(true$W, eStep$M); abline(a=0, b=1, h=0, v=0)
plot(true$Z, eStep$M%*%t(mStep$C)); abline(a=0, b=1, h=0, v=0)
lm(as.vector(mStep$C%*%t(mStep$C)) ~ -1 + as.vector(true$C%*%t(true$C)))$coef
