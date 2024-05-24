# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
# seed <- 1; set.seed(seed)
seed <- .Random.seed
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
cat(iter, ':', elboPath[iter], diff)
while((iter < iterMax) & (diff > tol)){
  iter <- iter+1
  mStepNew <- Mstep(data=data, mStep=mStep, eStep=eStep)
  eStepNew <- VEstep(data=data, mStep=mStepNew, eStep=eStep)
  elboPath[iter] <- ELBO(data=data, mStep=mStepNew, eStep=eStepNew)
  if(elboPath[iter] < elboPath[iter-1]){
    cat(iter, ':', 'E-1', ELBO(data=data, mStep=mStep, eStep=eStep), 
            'M', ELBO(data=data, mStep=mStepNew, eStep=eStep), 
            'E', ELBO(data=data, mStep=mStepNew, eStep=eStepNew), '\n')
  }
  diff <- max(max(abs(eStepNew$M - eStep$M)), max(abs(mStepNew$C - mStep$C)))
  mStep <- mStepNew; eStep <- eStepNew
  if(iter%%round(sqrt(iterMax))==0){
    plot(elboPath[1:iter], type='b', xlab='iter', ylim=quantile(elboPath[1:iter], probs=c(0.1, 1), na.rm=TRUE))
    cat(' /', iter, ':', elboPath[iter], diff)
  }
}
cat('\n')
elboPath <- elboPath[1:iter]
pseudoCovW <- t(eStep$M)%*%eStep$M/n + diag(colMeans(eStep$S))

par(mfrow=c(3, 3), pch=20, cex=0.75)
plot(elboPath, type='b', xlab='iter', ylim=quantile(elboPath, probs=c(0.1, 1)))
plot(true$gamma, mStep$gamma, ylim=range(c(mStep$gamma, oracle$mStep$gamma))); abline(a=0, b=1, h=0, v=0)
points(true$gamma, oracle$mStep$gamma, col=2)
points(true$gamma, init$mStep$gamma, col=4)
boxplot(eStep$xi ~ true$U)
plot(true$beta, mStep$beta); abline(a=0, b=1, h=0, v=0)
points(true$beta, oracle$mStep$beta, col=2)
points(true$beta, init$mStep$beta, col=4)
plot(true$C%*%t(true$C), mStep$C%*%t(mStep$C), ylim=range(cbind(mStep$C%*%t(mStep$C), oracle$mStep$C%*%t(oracle$mStep$C)))); abline(a=0, b=1, h=0, v=0)
points(true$C%*%t(true$C), mStep$C%*%pseudoCovW%*%t(mStep$C), col=3)
points(true$C%*%t(true$C), oracle$mStep$C%*%t(oracle$mStep$C), col=2);
points(true$C%*%t(true$C), init$mStep$C%*%t(init$mStep$C), col=4);
# plot(true$W, eStep$M); abline(a=0, b=1, h=0, v=0)
plot(true$Z, eStep$M%*%t(mStep$C)); abline(a=0, b=1, h=0, v=0)
plot(1+data$Y, 1+eStep$xi*NuMuA(data=data, mStep=mStep, eStep=eStep)$A, log='xy', xlab='Y', ylab='pred'); abline(a=0, b=1, h=0, v=0)

print(c(sum(diff(elboPath) < 0), 
        mean(diff(elboPath)[which(diff(elboPath) < 0)][-1]), 
        min(diff(elboPath))))

lm(as.vector(mStep$C%*%t(mStep$C)) ~ -1 + as.vector(true$C%*%t(true$C)))$coef
pseudoCovW
cov2cor(pseudoCovW)
