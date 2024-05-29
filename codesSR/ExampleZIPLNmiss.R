# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
source('Functions/FunctionsZIP.R')
source('Functions/FunctionsZIPLNmiss.R')
dataDir <- '../data/'
resDir <- '../results/'

# Parms
qList <- 1:10; qNb <- length(qList)

# Data
dataName <- 'Souchet'
# dataName <- 'Colvert'
X <- as.matrix(read.csv(paste0(dataDir, dataName, '_covariate.csv'), sep=',', header=TRUE))
Y <- as.matrix(read.csv(paste0(dataDir, dataName, '_count.csv'), sep=',', header=TRUE))
X <- X[, -1]; site <- Y[, 1]; Y <- Y[, -1]
print(c(dim(Y), dim(X)))
year <- as.numeric(gsub('X', '', colnames(Y)))
Omega <- 1*(!is.na(Y)); Y[which(Omega==0)] <- 0
data <- list(site=site, year=year, X=X, Y=Y, Omega=Omega, logFactY=lgamma(1+Y), 
             ij=cbind(rep(1:nrow(Y), ncol(Y)), rep(1:ncol(Y), each=nrow(Y))))
save(data, file=paste0(dataDir, dataName, '.Rdata'))

# Fit
for(qq in 1:qNb){
  q <- qList[qq]
  fitName <- paste0(dataName, '-ZiPLN-q', q)
  fitFile <- paste0(resDir, fitName, '.Rdata')
  if(!file.exists(fitFile)){
    print(fitFile)
    init <- InitZiPLN(data, q)
    vem <- VemZiPLN(data=data, init=init)
    save(init, vem, file=fitFile)
  }
}

