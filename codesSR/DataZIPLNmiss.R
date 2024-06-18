# Sim and fit ZI-PLN

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
library(bizicount); library(pscl)
source('Functions/FunctionsZIP.R')
source('Functions/FunctionsZIPLNmiss.R')
dataDir <- '../data/'
resDir <- '../results/'

# Data
dataName <- 'Souchet'
# dataName <- 'Colvert'
X <- as.matrix(read.csv(paste0(dataDir, dataName, '_covariate.csv'), sep=',', header=TRUE))
Y <- as.matrix(read.csv(paste0(dataDir, dataName, '_count.csv'), sep=',', header=TRUE))
X <- X[, -1]; site <- Y[, 1]; Y <- Y[, -1]
print(c(dim(Y), dim(X)))
year <- as.numeric(gsub('X', '', colnames(Y)))
Omega <- 1*(!is.na(Y)); Y[which(Omega==0)] <- 0

# Removing observations
obs <- 0.5; dataName <- paste0(dataName, '-obs', round(100*obs))
Omega <- matrix(rbinom(prod(dim(Y)), 1, obs), nrow(Y), ncol(Y)); Y[which(Omega==0)] <- 0

# Scaling covariates
Xmean <- colMeans(X); Xsd <- apply(X, 2, sd)
X <- as.matrix(cbind(rep(1, nrow(X)), scale(X)))
data <- list(site=site, year=year, X=X, Y=Y, Omega=Omega, 
             logFactY=lgamma(1+Y), Xmean=Xmean, Xsd=Xsd, 
             ij=cbind(rep(1:nrow(Y), ncol(Y)), rep(1:ncol(Y), each=nrow(Y))))
save(data, file=paste0(dataDir, dataName, '.Rdata'))

