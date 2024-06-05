# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
library(bizicount); library(pscl)
source('Functions/FunctionsZIP.R')
source('Functions/FunctionsZIPLNmiss.R')
source('Functions/FunctionsZIPLNmissVec.R')
dataDir <- '../data/'
resDir <- '../resultsSR/'

# Data
dataFull <- 'Colvert'; obs <- 0.5
load(paste0(dataDir, dataFull, '.Rdata'))
full <- data
dataName <- paste0(dataFull, '-obs', round(100*obs))
load(paste0(dataDir, dataName, '.Rdata'))
miss <- which(data$Omega==0)

# zip
zipFile <- paste0(resDir, dataName, '-zip.Rdata')
if(!file.exists(zipFile)){
  zip <- InitZiPLN(data=data, q=1)
  zip$mStep$C <- matrix(0, nrow(zip$mStep$C), ncol(zip$mStep$C))
  zip$eStep$M <- matrix(0, nrow(zip$eStep$M), ncol(zip$eStep$M))
  zip$eStep$S <- matrix(0, nrow(zip$eStep$S), ncol(zip$eStep$S))
  zipPred <- PredZiPLN(data=data, fit=zip)
  save(zip, zipPred, file=zipFile)
}else{load(zipFile)}

# Parms
qList <- 1:15;  qNb <- length(qList)

# Results
n <- nrow(data$Y); p <- ncol(data$Y); d <- ncol(data$X)
vemList <- list()
for(qq in 1:qNb){load(paste0(resDir, dataName, '-ZiPLN-q', qList[qq], '.Rdata'))
  vem$eStep$xi <- ComputeXi(data=data, mStep=vem$mStep, eStep=vem$eStep)
  vemList[[qq]] <- vem}
iter <- sapply(1:qNb, function(qq){vemList[[qq]]$iter})
diff <- sapply(1:qNb, function(qq){vemList[[qq]]$diff})
elbo <- sapply(1:qNb, function(qq){vemList[[qq]]$elbo})
penBic <- (2*d + choose(p, 2) - choose(qList-1, 2))*log(n)/2
bic <- elbo - penBic
ent <- 0.5*(sapply(1:qNb, function(qq){
  n*qList[qq]+log(2*pi*exp(1)) + sum(log(vemList[[qList[qq]]]$eStep$S))}))
ent <- ent - sapply(1:qNb, function(qq){
  xi <- vemList[[qList[qq]]]$eStep$xi; sum(xi*log(xi + (xi==0)) + (1 - xi)*log(1 - xi + (xi==1)))})
icl <- bic - ent
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
predList <- list()
for(qq in 1:qNb){load(paste0(resDir, dataName, '-ZiPLN-q', qList[qq], '-pred.Rdata'))
  predList[[qq]] <- pred}
par(mfrow=c(4, 4))
for(qq in 1:qNb){
  plot(1+predList[[qq]]$marg$esp[miss], 1+full$Y[miss], log='xy',
       xlab='predicted', ylab='observed', main=paste0('marginal q=', qList[qq]),
       col=1+(predList[[qq]]$marg$pi < 0.5)); abline(0, 1)
}
par(mfrow=c(4, 4))
for(qq in 1:qNb){
  plot(predList[[qq]]$cond$pi[miss], predList[[qq]]$cond$lambda[miss], log='y', 
       xlab='pi', ylab='lambda', main=paste0('conditional q=', qList[qq]), 
       col=1+(predList[[qq]]$cond$pi < 0.5))
}
par(mfrow=c(4, 4))
for(qq in 1:qNb){
  plot(1+predList[[qq]]$cond$esp[miss], 1+full$Y[miss], log='xy', 
       xlab='predicted', ylab='observed', main=paste0('conditional q=', qList[qq]), 
       col=1+(predList[[qq]]$cond$pi < 0.5)); abline(0, 1)
}
# par(mfrow=c(4, 4))
# for(qq in 1:qNb){
#   plot(1+predList[[qq]]$cond$esp, 1+full$Y[miss], log='xy', 
#        xlab='predicted', ylab='observed', main=paste0('conditional q=', qList[qq]), 
#        col=1+(predList[[qq]]$cond$pi < 0.5)); abline(0, 1)
#   points(1+predList[[qq]]$cond$lZipCI, 1+data$Y, col=4, pch='+'); abline(0, 1)
#   points(1+predList[[qq]]$cond$uZipCI, 1+data$Y, col=4, pch='+'); abline(0, 1)
# }
par(mfrow=c(4, 4))
for(qq in 1:qNb){
  boxplot(predList[[qq]]$cond$pi[miss] ~ (full$Y[miss]>0), 
       xlab='present', ylab='proba', main=paste0('conditional q=', qList[qq]))
}
piMat <- sapply(qList, function(q){predList[[q]]$cond$pi})
gammaMat <- sapply(qList, function(q){vemList[[q]]$mSte$gamma})
# plot(as.data.frame(gammaMat))

# Conditional moments
par(mfrow=c(1, 1))
qIcl <- which.max(icl)
vemList[[qIcl]]$eStep$MZ <- vemList[[qIcl]]$eStep$M%*%t(vemList[[qIcl]]$mStep$C)
vemList[[qIcl]]$eStep$SZ <- t(sapply(1:nrow(data$Y), function(i){diag(vemList[[qIcl]]$mStep$C%*%diag(vemList[[qIcl]]$eStep$S[i, ])%*%t(vemList[[qIcl]]$mStep$C))}))
plot(vemList[[qIcl]]$eStep$MZ, vemList[[qIcl]]$eStep$SZ, col=2-data$Omega)
boxplot(vemList[[qIcl]]$eStep$SZ ~ data$Omega)

# Compare with zip
par(mfrow=c(1, 1))
plot(1+full$Y[miss], 1+predList[[qIcl]]$cond$esp[miss], col=2, log='xy')
points(1+full$Y[miss], 1+zipPred$cond$esp[miss], col=8)
abline(0, 1, v=0, h=0)
plot(predList[[qIcl]]$marg$pi, zipPred$marg$pi); abline(0, 1, v=0, h=0)
plot(predList[[qIcl]]$marg$lambda, zipPred$marg$lambda, log='xy'); abline(0, 1, v=0, h=0)
