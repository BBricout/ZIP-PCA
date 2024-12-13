# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
seed <- 1; set.seed(seed)
# seed <- .Random.seed
source('Functions/FunctionsZIP.R')
source('Functions/FunctionsZIPLNmiss.R')
source('Functions/FunctionsUtils.R')
simDir <- '../simulSR/'
figDir <- '../plotsSR/'

# Parms
n <- 100; d <- 2; p <- 5; q <- 2; obs <- 0.95; coefC <- 0.1
# n <- 100; d <- 5; p <- 10; q <- 2; obs <- 0.6

# Simul
simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-coefC', coefC, '-seed', seed)
simNameFull <- paste0('ZiPLNsim', simParmsFull)
simFileFull <- paste0(simDir, simNameFull, '-noMiss.Rdata')
load(simFileFull)
simParms <- paste0(simParmsFull, '-obs', 100*obs)
simName <- paste0('ZiPLNsim', simParms)
simFile <- paste0(simDir, simName, '.Rdata')
if(!file.exists(simFile)){
  sim <- SimZiPLNmiss(sim=sim, obs=obs)
  save(sim, file=simFile)
}else{load(simFile)}
c(sum(diag(cov(matrix(data$X%*%true$mStep$beta, n, p)))), sum(diag(cov(true$latent$Z))))

# # First iterations
# ELBO(data=data, mStep=true$mStep, eStep=true$eStep)
# init <- InitZiPLN(data, q)
# ELBO(data=data, mStep=init$mStep, eStep=init$eStep)
# mStep <- mStep(data=data, mStep=init$mStep, eStep=init$eStep)
# ELBO(data=data, mStep=mStep, eStep=init$eStep)
# eStep <- VEstep(data=data, mStep=mStep, eStep=init$eStep)
# ELBO(data=data, mStep=mStep, eStep=eStep)

# Fit
fitName <- paste0('ZiPLNfit', simParms)
fitFile <- paste0(simDir, fitName, '.Rdata')
if(!file.exists(fitFile)){
  print(simName)
  init <- InitZiPLN(data=data, q=q)
  vem <- VemZiPLN(data, init=init)
  save(init, vem, file=fitFile)
}else{load(fitFile)}

# # Results
# oracle <- OracleZiPLN(data=data, latent=true$latent)
# # if(exportFig){png(paste0(figDir, fitName, '.png'))}
# par(mfrow=c(3, 3), pch=20, cex=0.75)
# plot(vem$elboPath, type='b', xlab='iter', ylim=quantile(vem$elboPath, probs=c(0.1, 1)))
# plot(true$mStep$gamma, vem$mStep$gamma, ylim=range(c(vem$mStep$gamma, oracle$mStep$gamma))); abline(a=0, b=1, h=0, v=0)
# points(true$mStep$gamma, oracle$mStep$gamma, col=2)
# points(true$mStep$gamma, init$mStep$gamma, col=8)
# boxplot(vem$eStep$xi ~ true$latent$U, col=2-data$Omega)
# hist(vem$eStep$xi, breaks=sqrt(n*p))
# plot(true$mStep$beta, vem$mStep$beta); abline(a=0, b=1, h=0, v=0)
# points(true$mStep$beta, oracle$mStep$beta, col=2)
# points(true$mStep$beta, init$mStep$beta, col=8)
# plot(true$mStep$C%*%t(true$mStep$C), vem$mStep$C%*%t(vem$mStep$C), ylim=range(cbind(vem$mStep$C%*%t(vem$mStep$C), oracle$mStep$C%*%t(oracle$mStep$C)))); abline(a=0, b=1, h=0, v=0)
# points(true$mStep$C%*%t(true$mStep$C), oracle$mStep$C%*%t(oracle$mStep$C), col=2);
# points(true$mStep$C%*%t(true$mStep$C), init$mStep$C%*%t(init$mStep$C), col=8);
# # plot(true$latent$W, vem$eStep$M); abline(a=0, b=1, h=0, v=0)
# plot(true$latent$Z, vem$eStep$M%*%t(vem$mStep$C)); abline(a=0, b=1, h=0, v=0)
# plot(1+true$latent$Yfull, 1+vem$pred$Yhat, log='xy', xlab='Y', ylab='pred', col=2-data$Omega); abline(a=0, b=1, h=0, v=0)
# points(1+true$latent$Yfull[which(data$Omega==0)], 1+vem$pred$Yhat[which(data$Omega==0)], col=2)
# # if(exportFig){dev.off()}

# Results
lm(as.vector(vem$mStep$C%*%t(vem$mStep$C)) ~ -1 + as.vector(true$mStep$C%*%t(true$mStep$C)))$coef
pseudoCovW <- t(vem$eStep$M)%*%vem$eStep$M/n + diag(colMeans(vem$eStep$S))
pseudoCovW
c(ELBO(data=data, mStep=true$mStep, eStep=true$eStep), 
  ELBO(data=data, mStep=init$mStep, eStep=true$eStep),
  ELBO(data=data, mStep=vem$mStep, eStep=vem$eStep))

# Jackknife
covName <- paste0('ZiPLNcov', simParms)
covFile <- paste0(simDir, covName, '.Rdata')
if(!file.exists(covFile)){
  jk <- JackknifeZiPLN(data=data, fit=vem)
  save(jk, file=covFile)
}else{load(covFile)}

par(mfrow=c(2, 2))
plot(vem$eStep$xi, jk$eStep$xi); abline(a=0, b=1, h=0, v=0)
plot(vem$eStep$M, jk$eStep$M); abline(a=0, b=1, h=0, v=0)
plot(vem$eStep$S, jk$eStep$S); abline(a=0, b=1, h=0, v=0)

vem$iter
summary(sapply(1:n, function(i){jk$fit_iList[[i]]$iter}))
vem$elbo
hist(sapply(1:n, function(i){diff(range(jk$fit_iList[[i]]$elboPath))}), main='', xlab='', ylab='')

rbind(true=true$mStep$gamma, vem=vem$mStep$gamma, sd=sqrt(diag(jk$cov$gamma)),
      stat=vem$mStep$gamma / sqrt(diag(jk$cov$gamma)), 
      statTrue=(vem$mStep$gamma - true$mStep$gamma) / sqrt(diag(jk$cov$gamma)))

rbind(true=true$mStep$beta, vem=vem$mStep$beta, sd=sqrt(diag(jk$cov$beta)),
      stat=vem$mStep$beta / sqrt(diag(jk$cov$beta)), 
      statTrue=(vem$mStep$beta - true$mStep$beta) / sqrt(diag(jk$cov$beta)))

# Predictions
