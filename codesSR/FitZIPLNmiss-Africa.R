# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

source('Functions/FunctionsUtils.R')
source('Functions/FunctionsZIP.R')
source('Functions/FunctionsZIPLNmiss.R')
splitMS <- TRUE

# Dirs 
dataDir <- '../dataSR/'; dataName <- 'Souchet-raw'
dataDir <- '../dataSR/'; dataName <- 'Souchet-miss50'
# dataDir <- '../dataSR/'; dataName <- 'Souchet-year2000-raw'
dataDir <- '../dataSR/'; dataName <- 'Souchet-year2000-miss50'
dataFile <- paste0(dataDir, dataName, '.Rdata')
load(dataFile)
n <- nrow(data$Y); d <- ncol(data$X); p <- ncol(data$Y); 

# Algo parms
orthC <- FALSE; q <- 2

# Fit
initFile <- paste0(dataDir, dataName, '-initZiPLNmiss.Rdata')
if(!file.exists(initFile)){
  print(dataName)
  init <- InitZiPLN(data=data, q=q, plot=TRUE)
  save(init, file=initFile)
}else{load(initFile)}
if(splitMS){
  fitFile <- paste0(dataDir, dataName, '-fitZiPLNmiss-splitMS.Rdata')
}else{
  fitFile <- paste0(dataDir, dataName, '-fitZiPLNmiss-splitMS', splitMS, '.Rdata')
}
print(fitFile)
if(!file.exists(fitFile)){
  print(dataName)
  vem <- try(VemZiPLN(data=data, init=init, orthC=orthC, splitMS=splitMS))
  # vem <- try(VemZiPLN(data=data, init=vem, orthC=orthC))
  save(init, vem, file=fitFile)
}else{load(fitFile)}

# Gradients
gradC <- ElboGradC(as.vector(vem$mStep$C), data=data, mStep=vem$mStep, eStep=vem$eStep)
gradM <- t(sapply(1:n, function(i){
  datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], 
                Omegai=data$Omega[i, ], logFactYi=data$logFactY[i, ])
  eStepi <- list(xii=vem$eStep$xi[i, ], mi=as.vector(vem$eStep$M[i, , drop=FALSE]), 
                 Si=as.vector(vem$eStep$S[i, , drop=FALSE]))
  ElboGradMi(mi=vem$eStep$M[i, ], datai=datai, mStep=vem$mStep, eStepi=eStepi)
}))
gradS <- t(sapply(1:n, function(i){
  datai <- list(Yi=data$Y[i, ], Xi=data$X[which(data$ij[, 1]==i), ], 
                Omegai=data$Omega[i, ], logFactYi=data$logFactY[i, ])
  eStepi <- list(xii=vem$eStep$xi[i, ], mi=as.vector(vem$eStep$M[i, , drop=FALSE]), 
                 Si=as.vector(vem$eStep$S[i, , drop=FALSE]))
  ElboGradSi(Si=vem$eStep$S[i, ], datai=datai, mStep=vem$mStep, eStepi=eStepi)
}))

par(mfrow=c(3, 3))
plot(vem$elboPath, type='b', xlab='iter', ylab='elbo', main='ZiPLN', 
     ylim=quantile(vem$elboPath, probs=c(0.1, 1), na.rm=TRUE))
hist(vem$eStep$M, breaks=sqrt(prod(dim(gradM))))
hist(vem$eStep$S, breaks=sqrt(prod(dim(gradS))))
hist(log10(abs(gradC)))
hist(log10(abs(gradM)), breaks=sqrt(prod(dim(gradM))))
hist(log10(abs(gradS)), breaks=sqrt(prod(dim(gradS))))

# Prediction
library(bizicount)
pred <- PredZiPLN(data=data, fit=vem)
plot(colMeans(data$Y, na.rm=TRUE), type='b')
points(colMeans(pred$marg$esp), type='b', col=4)
points(colMeans(pred$cond$esp), type='b', col=2)
points(colMeans(pred$nuMuA$A), type='b', col=3)

plot(data.frame(M=as.vector(vem$eStep$M), S=as.vector(vem$eStep$S), 
                log10gradM=as.vector(log10(abs(gradM))),
                log10gradS=as.vector(log10(abs(gradS)))))
