rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

#setwd("~/Documents/ZIP-PCA")

library(Rcpp)
library(PLNmodels)
library(missForest)

#---------- Functions SR
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIPLNmiss.R')
source('codesSR/Functions/FunctionsZIP.R')
#--- Codes BB
source("codesBB/FunctionsBB.R")
source("codesBB/UtilsBB.R")

###################################################"
#--- Simulate Données
#####################################################

# Parms: many small sims
n <- 1000; d <- 5; p <- 10; q <- 2
baseSimName <- 'ZiPLNsim'
seedList <- 1:10; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)

seedX <- 1; set.seed(seedX)
X0 <- matrix(rnorm(n*p*d), n*p, d); X0[, 1] <- 1; 

seed <- proc.time(); set.seed(seed)
sim <- SimZiPLN(n=n, p=p, d=d, q=q, X=X0)
data <- list(X=sim$X, Y=sim$Y,  ij=sim$ij, logFactY=lgamma(sim$Y+1))
true <- list(mStep=list(gamma=sim$gamma, beta=sim$beta, C=sim$C),
             eStep=list(M=sim$W, S=matrix(1e-4, n, q), xi=matrix(plogis(sim$X%*%sim$gamma), n, p)),
             latent=list(U=sim$U, W=sim$W, Z=sim$Z, Yall=sim$Yall))


tolXi <- 1e-04 ; tolS <- 1e-06

lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))

config <- PLNPCA_param()$config_optim
config$maxeval = 1000
# config$algorithm <- "MMA" # Par défaut dans config c'est "CCSAQ"
config$lower_bounds <- lb

##################### NON missing data
out.ref <- Miss.ZIPPCA(Y = data$Y, X = data$X, q, config = config)
#--------------- ELBO 
par(mfrow=c(1, 1))
plot(out.ref$elboPath[out.ref$elboPath > out.ref$elboPath[1]],type='l', main = "ELBO path",ylab='ELBO',xlab='iter')

############################# Missing data

#data$Y.na <- prodNA(data$Y, 0.01)

myrow <- sample(6:450,1) 
data$Y.na <- data$Y#prodNA(data$Y, 0.01)
#w <- which(data$Y[myrow,]==0)
data$Y.na[myrow,1] = NA
data$Y.na[(myrow-5):(myrow+5),]
data$Y[(myrow-5):(myrow+5),]


out.na <- Miss.ZIPPCA(Y = data$Y.na, X = data$X, q, config = config)

#----------------- remove the corresponding missing row
w <-1
databis <- list()
databis$Y <- data$Y[-myrow,]
databis$Y.na <- data$Y.na[-myrow,]
databis$logFactY <- data$logFactY[-myrow,]
pastis <- which(data$ij[,1] %in% myrow)
databis$X <- data$X[-pastis,]
databis$Omega <- data$Omega[-myrow,]
databis$R <- databis$Omega
config$lower_bounds <- c(rep(-Inf, 2*d + q*(p+n-length(myrow))), rep(tolS, (n-length(myrow)) * q))


out.na.databis <- Miss.ZIPPCA(Y = databis$Y.na, X = databis$X, q, config = config)
out.na.databis$elbo

plot(out.na$elboPath[out.na$elboPath > out.na$elboPath[1]],type='l', main = "ELBO path",ylab='ELBO',xlab='iter')
lines(out.ref$elboPath[out.ref$elboPath > out.ref$elboPath[1]],col='red')
lines(out.na.databis$elboPath[out.na.databis$elboPath > out.na.databis$elboPath[1]],col='green')


out.na.databis$elboPath





