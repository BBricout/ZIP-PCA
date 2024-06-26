# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
seed <- 1; set.seed(seed)
# seed <- .Random.seed
source('FunctionsZIP.R'); 
source('FunctionsZIPLN.R')
simDir <- '../simulSR/'
figDir <- '../plotsSR/'
exportFig <- TRUE

# Parms
n <- 500; d <- 20; p <- 30; q <- 2; obs <- 0.6

# Simul
simParms <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
simName <- paste0('ZiPLNsim', simParms)
simFile <- paste0(simDir, simName, '.Rdata')
if(!file.exists(simFile)){
  sim <- SimZiPLN(n=n, p=p, d=d, q=q)
  save(sim, file=simFile)
}else{load(simFile)}
data <- list(X=sim$X, Y=sim$Y, ij=sim$ij, logFactY=lgamma(sim$Y+1))
# sim <- SimZiPLN(n=n, p=p, d=d, q=q, obs=1)
# data <- list(X=sim$X, Y=sim$Y, Omega=sim$Omega, ij=sim$ij, logFactY=lgamma(sim$Y+1))

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

# Truth
true <- list(mstep=list(gamma=sim$gamma, beta=sim$beta, C=sim$C), 
             eStep=list(M=sim$W, S=matrix(1e-4, n, q), xi=matrix(plogis(sim$X%*%sim$gamma), n, p)), 
             latent=list(U=sim$U, W=sim$W, Z=sim$Z)             )

# Results
par(mfrow=c(3, 3), pch=20, cex=0.75)
if(exportFig){png(paste0(figDir, fitName, '.png'))}
plot(vem$elboPath, type='b', xlab='iter', ylim=quantile(vem$elboPath, probs=c(0.1, 1)))
plot(true$mStep$gamma, vem$mStep$gamma, ylim=range(c(vem$mStep$gamma, oracle$mStep$gamma))); abline(a=0, b=1, h=0, v=0)
points(true$mStep$gamma, oracle$mStep$gamma, col=2)
points(true$mStep$gamma, init$mStep$gamma, col=8)
boxplot(vem$eStep$xi ~ true$latent$U)
hist(vem$eStep$xi, breaks=sqrt(n*p))
plot(true$mStep$beta, vem$mStep$beta); abline(a=0, b=1, h=0, v=0)
points(true$mStep$beta, oracle$mStep$beta, col=2)
points(true$mStep$beta, init$mStep$beta, col=8)
plot(true$mStep$C%*%t(true$mStep$C), vem$mStep$C%*%t(vem$mStep$C), ylim=range(cbind(vem$mStep$C%*%t(vem$mStep$C), oracle$mStep$C%*%t(oracle$mStep$C)))); abline(a=0, b=1, h=0, v=0)
points(true$mStep$C%*%t(true$mStep$C), oracle$mStep$C%*%t(oracle$mStep$C), col=2);
points(true$mStep$C%*%t(true$mStep$C), init$mStep$C%*%t(init$mStep$C), col=8);
# plot(true$latent$W, vem$eStep$M); abline(a=0, b=1, h=0, v=0)
plot(true$latent$Z, vem$eStep$M%*%t(vem$mStep$C)); abline(a=0, b=1, h=0, v=0)
plot(1+data$Y, 1+vem$pred$Yhat, log='xy', xlab='Y', ylab='pred'); abline(a=0, b=1, h=0, v=0)
if(exportFig){dev.off()}

lm(as.vector(vem$mStep$C%*%t(vem$mStep$C)) ~ -1 + as.vector(true$mStep$C%*%t(true$mStep$C)))$coef
pseudoCovW <- t(vem$eStep$M)%*%vem$eStep$M/n + diag(colMeans(vem$eStep$S))
pseudoCovW
