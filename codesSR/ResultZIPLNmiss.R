# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

# seed <- .Random.seed
source('Functions/FunctionsZIPLNmiss.R')
library('bizicount')
simDir <- '../simulSR/'
figDir <- '../plotsSR/'
exportFig <- FALSE

# Parms
n <- 100; d <- 5; p <- 10; q <- 2
baseSimName <- 'ZiPLNsim'; baseFitName <- 'ZiPLNfit'; 
# baseSimName <- 'ZiPLNsim-sameX1'; baseFitName <- 'ZiPLNfit-sameX1';
seedList <- 1:10; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)

for(seed in seedList){ # seed <- 2
  for(oo in 1:obsNb){ # oo <- 3
    obs <- obsList[[oo]]
    # Data
    simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
    simNameFull <- paste0(baseSimName, simParmsFull)
    simFileFull <- paste0(simDir, simNameFull, '-noMiss.Rdata')
    simParms <- paste0(simParmsFull, '-obs', 100*obs)
    simName <- paste0(baseSimName, simParms)
    simFile <- paste0(simDir, simName, '.Rdata')
    fitName <- paste0(baseFitName, simParms)
    fitFile <- paste0(simDir, fitName, '.Rdata')
    if(file.exists(fitFile)){
      print(fitName)
      load(simFileFull)
      true$latent$Yfull <- data$Y
      load(simFile)
      load(fitFile)
      # Oracle & preds
      oracle <- OracleZiPLN(data=data, latent=true$latent)
      pred <- PredZiPLN(data=data, fit=vem)
      # Truth
      elboTrue <- ELBO(data=data, mStep=true$mstep, eStep=true$eStep)
      elboInit <- ELBO(data=data, mStep=init$mStep, eStep=init$eStep)
      # Results
      if(exportFig){png(paste0(figDir, fitName, '.png'))}
      par(mfrow=c(3, 3), pch=20, cex=0.75)
      plot(vem$elboPath, type='b', xlab='iter', main=simParms, ylim=quantile(vem$elboPath, probs=c(0.1, 1)))
      abline(h=elboTrue, col=2, lty=2, lwd=2)
      text(0.6*vem$iter, mean(quantile(vem$elboPath, probs=c(0.1, 1))), label=round(elboTrue))
      text(0.8*vem$iter, mean(quantile(vem$elboPath, probs=c(0.1, 1))), label=round(elboInit))
      plot(true$mstep$gamma, init$mStep$gamma, col=4, pch=21, ylim=range(c(vem$mStep$gamma, oracle$mStep$gamma))); abline(a=0, b=1, h=0, v=0)
      points(true$mstep$gamma, oracle$mStep$gamma, col=2, pch=21)
      points(true$mstep$gamma, vem$mStep$gamma)
      boxplot(vem$eStep$xi ~ true$latent$U)
      hist(vem$eStep$xi, breaks=sqrt(n*p))
      plot(true$mstep$beta, init$mStep$beta, col=4, pch=21); abline(a=0, b=1, h=0, v=0)
      points(true$mstep$beta, oracle$mStep$beta, col=2, pch=21)
      points(true$mstep$beta, vem$mStep$beta)
      plot(true$mstep$C%*%t(true$mstep$C), init$mStep$C%*%t(init$mStep$C), col=4, pch=21, ylim=range(cbind(vem$mStep$C%*%t(vem$mStep$C), oracle$mStep$C%*%t(oracle$mStep$C)))); abline(a=0, b=1, h=0, v=0)
      points(true$mstep$C%*%t(true$mstep$C), oracle$mStep$C%*%t(oracle$mStep$C), col=2, pch=21);
      points(true$mstep$C%*%t(true$mstep$C), vem$mStep$C%*%t(vem$mStep$C));
      # plot(true$W, vem$eStep$M); abline(a=0, b=1, h=0, v=0)
      plot(true$latent$Z, vem$eStep$M%*%t(vem$mStep$C)); abline(a=0, b=1, h=0, v=0)
      plot(1+true$latent$Yfull, 1+vem$pred$Yhat, log='xy', xlab='Y', ylab='pred', col=2-data$Omega); abline(a=0, b=1, h=0, v=0)
      points(1+true$latent$Yfull[which(data$Omega==0)], 1+vem$pred$Yhat[which(data$Omega==0)], col=2)
      if(exportFig){dev.off()}
      
      lm(as.vector(vem$mStep$C%*%t(vem$mStep$C)) ~ -1 + as.vector(true$mstep$C%*%t(true$mstep$C)))$coef
      pseudoCovW <- t(vem$eStep$M)%*%vem$eStep$M/n + diag(colMeans(vem$eStep$S))
      pseudoCovW
      
      miss <- which(data$Omega==0); missNb <- length(miss)
      # Presence
      margPres <- matrix(1*(pred$marg$pi>.5), n, p); 
      condPres <- matrix(1*(pred$cond$pi>.5), n, p)
      print(c(sum(diag(table(true$latent$U[miss], margPres[miss])))/missNb,
              sum(diag(table(true$latent$U[miss], condPres[miss])))/missNb))
      # Abundance: CI cover
      print(c(mean(((true$latent$Yfull - pred$marg$lZipCI)*(pred$marg$uZipCI - true$latent$Yfull) > 0)[miss]), 
              mean(((true$latent$Yfull - pred$cond$lZipCI)*(pred$cond$uZipCI - true$latent$Yfull) > 0)[miss])))
      # Abundance: CI cover conditional on presence
      missMargPres <- which((data$Omega==0) & (margPres==1))
      missCondPres <- which((data$Omega==0) & (condPres==1))
      print(c(mean(((true$latent$Yfull - pred$marg$lZipCI)*(pred$marg$uZipCI - true$latent$Yfull) > 0)[missMargPres]), 
              mean(((true$latent$Yfull - pred$cond$lZipCI)*(pred$cond$uZipCI - true$latent$Yfull) > 0)[missCondPres])))
    }
  }
}
