# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
library(bizicount); library(pscl)
source('Functions/FunctionsZIP.R')
source('Functions/FunctionsUtils.R')
source('Functions/FunctionsZIPLNpost.R')
source('Functions/FunctionsZIPLNmiss.R')
dataDir <- '../data/'
resDir <- '../resultsSR/'

# Data
dataName <- 'Souchet'
dataName <- 'Colvert'
obs <- 0.5; dataName <- paste0(dataName, '-obs', round(100*obs))
load(paste0(dataDir, dataName, '.Rdata'))

# Parms
qList <- 1:ncol(data$Y); qNb <- length(qList)
qList <- 1:15;  qNb <- length(qList)

# Fit
for(qq in 1:qNb){
  q <- qList[qq]
  fitName <- paste0(dataName, '-ZiPLN-q', q)
  fitFile <- paste0(resDir, fitName, '.Rdata')
  if(!file.exists(fitFile)){
    print(fitFile)
    init <- InitZiPLN(data, q)
    vem <- VemZiPLN(data=data, init=init, iterMax=5e3)
    save(init, vem, file=fitFile)
  }
}

# Results
n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X)
vemList <- list()
for(qq in 1:qNb){load(paste0(resDir, dataName, '-ZiPLN-q', qList[qq], '.Rdata'))
  vem$eStep$xi <- ComputeXi(data=data, mStep=vem$mStep, eStep=vem$eStep)
  vemList[[qq]] <- vem}
critList <- lapply(vemList, function(vem){Criteria(data=data, fit=vem)})
elbo <- unlist(lapply(critList, function(crit){crit$elbo}))
bic <- unlist(lapply(critList, function(crit){crit$bic}))
icl <- unlist(lapply(critList, function(crit){crit$icl}))
plot(qList, elbo, type='b', ylim=quantile(c(elbo, bic, icl), prob=c(.1, 1)))
points(qList, bic, type='b', col=2)
points(qList, icl, type='b', col=4)
abline(v=qList[which.max(bic)], col=2, lty=2)
abline(v=qList[which.max(icl)], col=4, lty=2)

# # Results
# qq <- qList[which.max(icl)]
# jkName <- paste0(dataName, '-ZiPLN-q', qList[qq], '-jackknife')
# jkFile <- paste0(resDir, jkName, '.Rdata')
# if(!file.exists(jkFile)){
#   print(jkFile)
#   jk <- JackknifeZiPLN(data=data, fit=vemList[[qq]])
#   save(jk, file=jkFile)
# }

# Preds
for(qq in 1:qNb){
  q <- qList[qq]
  predName <- paste0(dataName, '-ZiPLN-q', q, '-pred')
  predFile <- paste0(resDir, predName, '.Rdata')
  if(!file.exists(predFile)){
    print(predFile)
    pred <- PredZiPLN(data=data, fit=vemList[[qq]])
    save(pred, file=predFile)
  }
}
predList <- list()
for(qq in 1:qNb){load(paste0(resDir, dataName, '-ZiPLN-q', qList[qq], '-pred.Rdata'))
  predList[[qq]] <- pred}
# par(mfrow=c(4, 3))
# for(qq in 1:qNb){
#   plot(1+predList[[qq]]$marg$esp, 1+data$Y, log='xy', 
#        xlab='predicted', ylab='observed', main=paste0('marginal q=', qList[qq]), 
#        col=1+(predList[[qq]]$marg$pi < 0.5)); abline(0, 1)
# }
par(mfrow=c(4, 4))
for(qq in 1:qNb){
  plot(predList[[qq]]$cond$pi, predList[[qq]]$cond$lambda, log='y', 
       xlab='pi', ylab='lambda', main=paste0('conditional q=', qList[qq]), 
       col=1+(predList[[qq]]$cond$pi < 0.5))
}
par(mfrow=c(4, 4))
for(qq in 1:qNb){
  plot(1+predList[[qq]]$cond$esp, 1+data$Y, log='xy', 
       xlab='predicted', ylab='observed', main=paste0('conditional q=', qList[qq]), 
       col=1+(predList[[qq]]$cond$pi < 0.5)); abline(0, 1)
}
par(mfrow=c(4, 4))
for(qq in 1:qNb){
  plot(1+predList[[qq]]$cond$esp, 1+data$Y, log='xy', 
       xlab='predicted', ylab='observed', main=paste0('conditional q=', qList[qq]), 
       col=1+(predList[[qq]]$cond$pi < 0.5)); abline(0, 1)
  points(1+predList[[qq]]$cond$esp, 1+predList[[qq]]$cond$lZipCI, col=4, pch='+'); abline(0, 1)
  points(1+predList[[qq]]$cond$esp, 1+predList[[qq]]$cond$uZipCI, col=4, pch='+'); abline(0, 1)
}

