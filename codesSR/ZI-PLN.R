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

# Fit
oracle <- OracleZiPLN(sim)
init <- InitZiPLN(data); mStep <- init$mStep; eStep <- init$eStep
vem <- VemZiPLN(data, init=init)

# Results
par(mfrow=c(3, 3), pch=20, cex=0.75)
plot(elboPath, type='b', xlab='iter', ylim=quantile(elboPath, probs=c(0.1, 1)))
plot(true$gamma, vem$mStep$gamma, ylim=range(c(vem$mStep$gamma, oracle$mStep$gamma))); abline(a=0, b=1, h=0, v=0)
points(true$gamma, oracle$mStep$gamma, col=2)
points(true$gamma, init$mStep$gamma, col=4)
boxplot(vem$eStep$xi ~ true$U)
plot(true$beta, vem$mStep$beta); abline(a=0, b=1, h=0, v=0)
points(true$beta, oracle$mStep$beta, col=2)
points(true$beta, init$mStep$beta, col=8)
plot(true$C%*%t(true$C), vem$mStep$C%*%t(vem$mStep$C), ylim=range(cbind(vem$mStep$C%*%t(vem$mStep$C), oracle$mStep$C%*%t(oracle$mStep$C)))); abline(a=0, b=1, h=0, v=0)
points(true$C%*%t(true$C), oracle$mStep$C%*%t(oracle$mStep$C), col=2);
points(true$C%*%t(true$C), init$mStep$C%*%t(init$mStep$C), col=8);
# plot(true$W, vem$eStep$M); abline(a=0, b=1, h=0, v=0)
plot(true$Z, vem$eStep$M%*%t(vem$mStep$C)); abline(a=0, b=1, h=0, v=0)
plot(1+data$Y, 1+vem$pred$Yhat, log='xy', xlab='Y', ylab='pred'); abline(a=0, b=1, h=0, v=0)

print(c(sum(diff(elboPath) < 0), 
        mean(diff(elboPath)[which(diff(elboPath) < 0)][-1]), 
        min(diff(elboPath))))

lm(as.vector(mStep$C%*%t(mStep$C)) ~ -1 + as.vector(true$C%*%t(true$C)))$coef
pseudoCovW <- t(eStep$M)%*%eStep$M/n + diag(colMeans(eStep$S))
pseudoCovW
