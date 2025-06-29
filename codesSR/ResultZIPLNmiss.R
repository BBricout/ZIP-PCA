# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

# seed <- .Random.seed
source('Functions/FunctionsZIPLNmiss.R')
library('bizicount')
simDir <- '../simulSR/'
figDir <- '../plotsSR/'
exportFig <- FALSE


# Parms
n <- 1000; d <- 2; p <- 5; q <- 2; coefC <- 1
vec <- FALSE
baseSimName <- 'ZiPLNsim'; baseFitName <- 'ZiPLNfit'; basePlotName <- 'ZiPLNplot'
# baseSimName <- 'ZiPLNsim-sameX1'; baseFitName <- 'ZiPLNfit-sameX1';
# seedList <- 1:10; seedNb <- length(seedList)
# obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)
seedList <- 1:100; seedNb <- length(seedList)
obsList <- c(1); obsNb <- length(obsList)
parmsGene <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-coefC', coefC)

for(oo in 1:obsNb){ # oo <- 1
  obs <- obsList[[oo]]
  thetaTrue <- thetaHat <- matrix(NA, seedNb, 2*d+p*q)
  colnames(thetaTrue) <- colnames(thetaHat) <- c(paste0('gamma', 1:d), paste0('beta', 1:d), paste0('C', 1:(p*q)))
  iter <- rep(NA, seedNb)
  for(seed in seedList){ # seed <- 1
    # Data
    parmsName <- paste0(parmsGene, '-seed', seed)
    simNameFull <- paste0(baseSimName, parmsName)
    simFileFull <- paste0(simDir, simNameFull, '-noMiss.Rdata')
    simParms <- paste0(parmsName, '-obs', 100*obs)
    simName <- paste0(baseSimName, simParms)
    simFile <- paste0(simDir, simName, '.Rdata')
    fitName <- paste0(baseFitName, simParms)
    fitFile <- paste0(simDir, fitName, '.Rdata')
    if(vec){fitFile <- paste0(simDir, fitName, '-vec.Rdata')}
    if(file.exists(fitFile)){
      print(fitName)
      load(simFileFull)
      true$latent$Yfull <- data$Y
      # load(simFile)
      thetaTrue[seed, ] <- c(true$mStep$gamma, true$mStep$beta, as.vector(true$mStep$C))
      load(fitFile)
      if(length(vem) > 1){
        thetaHat[seed, ] <- c(vem$mStep$gamma, vem$mStep$beta, as.vector(vem$mStep$C))
        iter[seed] <- vem$iter
        cat(thetaTrue[seed, 1:(2*d)], thetaHat[seed, 1:(2*d)], '\n')
      }
      # # Oracle & preds
      # oracle <- OracleZiPLN(data=data, latent=true$latent)
      # pred <- PredZiPLN(data=data, fit=vem)
      # # Truth
      # elboTrue <- ELBO(data=data, mStep=true$mStep, eStep=true$eStep)
      # elboInit <- ELBO(data=data, mStep=init$mStep, eStep=init$eStep)
      
      # Results
      # if(exportFig){png(paste0(figDir, fitName, '.png'))}
      # par(mfrow=c(3, 3), pch=20, cex=0.75)
      # plot(vem$elboPath, type='b', xlab='iter', main=simParms, ylim=quantile(vem$elboPath, probs=c(0.1, 1)))
      # abline(h=elboTrue, col=2, lty=2, lwd=2)
      # text(0.6*vem$iter, mean(quantile(vem$elboPath, probs=c(0.1, 1))), label=round(elboTrue))
      # text(0.8*vem$iter, mean(quantile(vem$elboPath, probs=c(0.1, 1))), label=round(elboInit))
      # plot(true$mStep$gamma, init$mStep$gamma, col=4, pch=21, ylim=range(c(vem$mStep$gamma, oracle$mStep$gamma))); abline(a=0, b=1, h=0, v=0)
      # points(true$mStep$gamma, oracle$mStep$gamma, col=2, pch=21)
      # points(true$mStep$gamma, vem$mStep$gamma)
      # boxplot(vem$eStep$xi ~ true$latent$U)
      # hist(vem$eStep$xi, breaks=sqrt(n*p))
      # plot(true$mStep$beta, init$mStep$beta, col=4, pch=21); abline(a=0, b=1, h=0, v=0)
      # points(true$mStep$beta, oracle$mStep$beta, col=2, pch=21)
      # points(true$mStep$beta, vem$mStep$beta)
      # plot(true$mStep$C%*%t(true$mStep$C), init$mStep$C%*%t(init$mStep$C), col=4, pch=21, ylim=range(cbind(vem$mStep$C%*%t(vem$mStep$C), oracle$mStep$C%*%t(oracle$mStep$C)))); abline(a=0, b=1, h=0, v=0)
      # points(true$mStep$C%*%t(true$mStep$C), oracle$mStep$C%*%t(oracle$mStep$C), col=2, pch=21);
      # points(true$mStep$C%*%t(true$mStep$C), vem$mStep$C%*%t(vem$mStep$C));
      # # plot(true$W, vem$eStep$M); abline(a=0, b=1, h=0, v=0)
      # plot(true$latent$Z, vem$eStep$M%*%t(vem$mStep$C)); abline(a=0, b=1, h=0, v=0)
      # plot(1+true$latent$Yfull, 1+vem$pred$Yhat, log='xy', xlab='Y', ylab='pred', col=2-data$Omega); abline(a=0, b=1, h=0, v=0)
      # points(1+true$latent$Yfull[which(data$Omega==0)], 1+vem$pred$Yhat[which(data$Omega==0)], col=2)
      # if(exportFig){dev.off()}
      
      # lm(as.vector(vem$mStep$C%*%t(vem$mStep$C)) ~ -1 + as.vector(true$mStep$C%*%t(true$mStep$C)))$coef
      # pseudoCovW <- t(vem$eStep$M)%*%vem$eStep$M/n + diag(colMeans(vem$eStep$S))
      # pseudoCovW
      
      # miss <- which(data$Omega==0); missNb <- length(miss)
      # # Presence
      # margPres <- matrix(1*(pred$marg$pi>.5), n, p); 
      # condPres <- matrix(1*(pred$cond$pi>.5), n, p)
      # print(c(sum(diag(table(true$latent$U[miss], margPres[miss])))/missNb,
      #         sum(diag(table(true$latent$U[miss], condPres[miss])))/missNb))
      # # Abundance: CI cover
      # print(c(mean(((true$latent$Yfull - pred$marg$lZipCI)*(pred$marg$uZipCI - true$latent$Yfull) > 0)[miss]), 
      #         mean(((true$latent$Yfull - pred$cond$lZipCI)*(pred$cond$uZipCI - true$latent$Yfull) > 0)[miss])))
      # # Abundance: CI cover conditional on presence
      # missMargPres <- which((data$Omega==0) & (margPres==1))
      # missCondPres <- which((data$Omega==0) & (condPres==1))
      # print(c(mean(((true$latent$Yfull - pred$marg$lZipCI)*(pred$marg$uZipCI - true$latent$Yfull) > 0)[missMargPres]), 
      #         mean(((true$latent$Yfull - pred$cond$lZipCI)*(pred$cond$uZipCI - true$latent$Yfull) > 0)[missCondPres])))
    }
  }
  if(exportFig){png(paste0(figDir, basePlotName, parmsGene, '-obs', 100*obs, '.png'))}
  par(mfrow=c(2, 2))
  thetaSd <- apply(thetaHat, 2, sd, na.rm=TRUE)
  thetaStat <- (thetaHat - thetaTrue)/(rep(1, seedNb)%o%thetaSd)
  for(h in 1:4){qqnorm(thetaStat[, h], main=paste(coefC, obs, colnames(thetaTrue)[h]), na.rm=TRUE, pch=20); 
    abline(a=0, b=1, h=0, v=0)}
  if(exportFig){dev.off()}
  res <- rbind(colMeans((thetaHat - thetaTrue)[, 1+c(0, d)], na.rm=TRUE), 
               colMeans(thetaStat[, 1+c(0, d)], na.rm=TRUE), 
               thetaSd[1+c(0, d)])
  rownames(res) <- c('bias', 'stat', 'sd')
  print(res)
  print(summary(iter)); boxplot(iter)
}

# # Check pred
# plot(data$X%*%true$mStep$beta + as.vector(true$latent$W%*%t(true$mStep$C)), 
#      data$X%*%vem$mStep$beta + as.vector(vem$eStep$M%*%t(vem$mStep$C)))
# abline(0, 1, h=0, v=0)
# plot(as.vector(true$latent$W%*%t(true$mStep$C)), 
#      as.vector(vem$eStep$M%*%t(vem$mStep$C)))
# abline(0, 1, h=0, v=0)
