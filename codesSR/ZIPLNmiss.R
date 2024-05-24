# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
seed <- 4; set.seed(seed)
# seed <- .Random.seed
source('FunctionsZIP.R')
source('FunctionsZIPLNmiss.R')
simDir <- '../simulSR/'
figDir <- '../plotsSR/'
exportFig <- TRUE

# Parms
n <- 100; d <- 5; p <- 10; q <- 2; obs <- 0.9

# Simul
simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
simNameFull <- paste0('ZiPLNsim', simParmsFull)
simFileFull <- paste0(simDir, simNameFull, '.Rdata')
load(simFileFull)
simParms <- paste0(simParmsFull, '-obs', 100*obs)
simName <- paste0('ZiPLNsim', simParms)
simFile <- paste0(simDir, simName, '.Rdata')
if(!file.exists(simFile)){
  sim <- SimZiPLNmiss(sim=sim, obs=obs)
  save(sim, file=simFile)
}else{load(simFile)}
data <- list(X=sim$X, Y=sim$Y, Omega=sim$Omega, ij=sim$ij, logFactY=lgamma(sim$Y+1))
true <- list(gamma=sim$gamma, beta=sim$beta, C=sim$C, 
             U=sim$U, W=sim$W, Z=sim$Z,
             M=sim$W, S=matrix(1e-4, n, q), xi=matrix(plogis(sim$X%*%sim$gamma), n, p))
c(sum(diag(cov(matrix(sim$X%*%sim$beta, n, p)))), sum(diag(cov(sim$Z))))

# # First iterations
# init <- InitZiPLN(data)
# ELBO(data=data, mStep=init$mStep, eStep=init$eStep)
# mStep <- Mstep(data=data, mStep=init$mStep, eStep=init$eStep)
# ELBO(data=data, mStep=mStep, eStep=init$eStep)
# eStep <- VEstep(data=data, mStep=mStep, eStep=init$eStep)
# ELBO(data=data, mStep=mStep, eStep=eStep)

# Fit
fitName <- paste0('ZiPLNfit', simParms)
fitFile <- paste0(simDir, fitName, '.Rdata')
if(!file.exists(fitFile)){
  print(simName)
  oracle <- OracleZiPLN(sim)
  init <- InitZiPLN(data)
  vem <- VemZiPLN(data, init=init)
  save(oracle, init, vem, file=fitFile)
}else{load(fitFile)}

# Results
if(exportFig){png(paste0(figDir, fitName, '.png'))}
par(mfrow=c(3, 3), pch=20, cex=0.75)
plot(vem$elboPath, type='b', xlab='iter', ylim=quantile(vem$elboPath, probs=c(0.1, 1)))
plot(true$gamma, vem$mStep$gamma, ylim=range(c(vem$mStep$gamma, oracle$mStep$gamma))); abline(a=0, b=1, h=0, v=0)
points(true$gamma, oracle$mStep$gamma, col=2)
points(true$gamma, init$mStep$gamma, col=8)
boxplot(vem$eStep$xi ~ true$U, col=2-data$Omega)
hist(vem$eStep$xi, breaks=sqrt(n*p))
plot(true$beta, vem$mStep$beta); abline(a=0, b=1, h=0, v=0)
points(true$beta, oracle$mStep$beta, col=2)
points(true$beta, init$mStep$beta, col=8)
plot(true$C%*%t(true$C), vem$mStep$C%*%t(vem$mStep$C), ylim=range(cbind(vem$mStep$C%*%t(vem$mStep$C), oracle$mStep$C%*%t(oracle$mStep$C)))); abline(a=0, b=1, h=0, v=0)
points(true$C%*%t(true$C), oracle$mStep$C%*%t(oracle$mStep$C), col=2);
points(true$C%*%t(true$C), init$mStep$C%*%t(init$mStep$C), col=8);
# plot(true$W, vem$eStep$M); abline(a=0, b=1, h=0, v=0)
plot(true$Z, vem$eStep$M%*%t(vem$mStep$C)); abline(a=0, b=1, h=0, v=0)
plot(1+sim$Yfull, 1+vem$pred$Yhat, log='xy', xlab='Y', ylab='pred', col=2-data$Omega); abline(a=0, b=1, h=0, v=0)
points(1+sim$Yfull[which(data$Omega==0)], 1+vem$pred$Yhat[which(data$Omega==0)], col=2)
if(exportFig){dev.off()}

lm(as.vector(vem$mStep$C%*%t(vem$mStep$C)) ~ -1 + as.vector(true$C%*%t(true$C)))$coef
pseudoCovW <- t(vem$eStep$M)%*%vem$eStep$M/n + diag(colMeans(vem$eStep$S))
pseudoCovW
