# Reshape data for fitting ZiPLNmiss

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

library(lori)

# Dirs 
dataDir <- '../../Afrique/'
dataDirSR <- '../dataSR/'

# Counts
dataName <- 'Souchet'; Y <- read.table(paste0(dataDir, 'souchet.csv'), sep=';', header=TRUE)[, -1]
# dataName <- 'SarcelleLori'; Y <- read.table(paste0(dataDir, 'sarcelle_hiver_filtrageLori.csv'), sep=';', header=TRUE)[, -1]
head(Y)
Y <- as.matrix(Y); dim(Y); n <- nrow(Y); p <- ncol(Y)
years <- as.numeric(gsub('X', '', colnames(Y)))
yearThresh <- min(years)
yearThresh <- 2000
hist(colMeans(is.na(Y)), breaks=sqrt(n))
missThresh <- 50

# Covariates
covSite <- read.table(paste0(dataDir, 'cov_sites.csv'), sep=';', header=TRUE, dec=',')[, -1]
dim(covSite)
covYear <- read.table(paste0(dataDir, 'cov_years.csv'), sep=';', header=TRUE, dec=',')[, -1]
dim(covYear)
covSiteYear <- read.table(paste0(dataDir, 'cov_sites_years.csv'), sep=';', header=TRUE, dec=',')
dim(covSiteYear)
X <- covmat(n=n, p=p, R=covSite, C=covYear, E=covSiteYear)
X <- scale(X)

# Remove years
if(yearThresh > min(years)){
  jSel <- which(years >= yearThresh)
  ij <- cbind(rep(1:n, p), rep(1:p, each=n))
  ijSel <- which(ij[, 2]%in%jSel)
  X <- X[ijSel, ]; Y <- Y[, jSel]; dim(Y); n <- nrow(Y); p <- ncol(Y)
  X <- X[, which(apply(X, 2, sd) > 0)]
  dataName <- paste0(dataName, '-year', yearThresh)
}

# Remove poor sites
iSel <- which(rowMeans(is.na(Y)) <= missThresh/100)
ij <- cbind(rep(1:n, p), rep(1:p, each=n))
ijSel <- which(ij[, 1]%in%iSel)
X <- X[ijSel, ]; Y <- Y[iSel, ]; dim(Y); n <- nrow(Y); p <- ncol(Y)
X <- X[, which(apply(X, 2, sd) > 0)]
if(missThresh < 100){
  dataName <- paste0(dataName, '-miss', missThresh)
  }else{dataName <- paste0(dataName, '-raw')}

# Export
X <- as.matrix(cbind(rep(1, n*p), X))
Omega <- 1*(!is.na(Y))
Y[is.na(Y)] <- 0
ij <- cbind(rep(1:n, p), rep(1:p, each=n))
logFactY <- lgamma(1+Y)
logFactY[is.na(logFactY)] <- 0
data <- list(X=X, Y=Y, Omega=Omega, ij=ij, logFactY=logFactY)
save(data, file=paste0(dataDirSR, dataName, '.Rdata'))
