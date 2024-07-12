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


#--- Tests elbo

tolXi <- 1e-04 ; tolS <- 1e-06

# data.na <- data
# 
# # data.na$Y.na <- prodNA(data$Y, 0.01)
# 
# # data.na$Omega <- data.na$R <- ifelse(is.na(data.na$Y.na), 0, 1)
# data.na$Omega <- data.na$R <- data$R
# 
# data.na$Omega[1,1] <- data.na$R[1,1] <- 0
# 
# init <- InitZiPLN(data, q)
# 
# params <- list(B = as.matrix(init$mStep$beta), 
#                D = as.matrix(init$mStep$gamma), 
#                C = as.matrix(init$mStep$C), 
#                M = as.matrix(init$eStep$M), 
#                S = as.matrix(init$eStep$S))
# 
# # Tests xi
# 
# xiSR <- ComputeXi(data, init$mStep, init$eStep)
# xiBB <- Elbo_grad_Rcpp(data, params, tolXi)$xi
# 
# plot(xiSR, xiBB) ; abline(0,1)
# 
# #--- Tests Elbo
# 
# # Données complètes
# Belbo <- Elbo_grad_Rcpp(data, params, tolXi)
# Selbo <- ELBO(data, init$mStep, init$eStep)
# 
# elbo.test <- c(Belbo$elbo1, Belbo$elbo2, Belbo$elbo3, Belbo$elbo4, Belbo$elbo5, Belbo$objective)
# Selbo
# elbo.test
# 
# # Données manquantes 
# 
# init.na <- InitZiPLN(data.na, q, tolXi)
# params.na <- list(B = as.matrix(init.na$mStep$beta), 
#                D = as.matrix(init.na$mStep$gamma), 
#                C = as.matrix(init.na$mStep$C), 
#                M = as.matrix(init.na$eStep$M), 
#                S = as.matrix(init.na$eStep$S))
# 
# Belbo.na <- Elbo_grad_Rcpp(data.na, params.na, tolXi)
# Selbo.na <- Selbo <- ELBO(data.na, init.na$mStep, init.na$eStep)
# 
# elbo.test.na <- c(Belbo.na$elbo1, Belbo.na$elbo2, Belbo.na$elbo3, Belbo.na$elbo4, Belbo.na$elbo5, Belbo.na$objective)
# Selbo.na
# elbo.test.na
# 
# 
# 
# elbo.test
# elbo.test.na
# 
# plot(Belbo$xi, Belbo.na$xi) ; abline(0,1)



#--- Log S


lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))


config <- PLNPCA_param()$config_optim
config$maxeval = 1000
# config$algorithm <- "MMA" # Par défaut dans config c'est "CCSAQ"
config$lower_bounds <- lb
config$maxeval <- 200



out.ref <- Miss.ZIPPCA(Y = data$Y, X = data$X, q, config = config)
out.log <- Miss.ZIPPCA.logS(Y = data$Y, X = data$X, q, config = config)


# params.test <- Init_ZIP(prodNA(data$Y, 0.01), X = data$X, q = q)
# params.test$S <- matrix(0.1, n, q)

Y.test <- prodNA(data$Y, 0.1)
# Miss.ZIPPCA.logS(Y = Y.test, X = data$X, q, config = config)$elbo
# Miss.ZIPPCA.logS(Y = Y.test, X = data$X, q, params = params.test, config = config)$elbo


out.na <- Miss.ZIPPCA(Y = Y.test, X = data$X, q, config = config)
out.log.na <- Miss.ZIPPCA.logS(Y = Y.test, X = data$X, q, config = config)


# out <- out.ref


# # Estimations
# B.hat <- out$mStep$beta
# D.hat <- out$mStep$gamma
# C.hat <- out$mStep$C
# M.hat <- out$eStep$M
# S.hat <- out$eStep$S
# 
# XB <- VectorToMatrix(data$X %*% true$mStep$beta, n, p)
# XB.hat <- VectorToMatrix(data$X %*% B.hat, n, p)
# XB.init <- VectorToMatrix(data$X %*% out$params.init$B, n, p)
# 
# 
# XD <- VectorToMatrix(data$X %*% true$mStep$gamma, n, p)
# XD.hat <- VectorToMatrix(data$X %*% D.hat, n, p)
# XD.init <- VectorToMatrix(data$X %*% out$params.init$D, n, p)
# 
# 
# pred.init <- pred <- exp(XB.init) #  + params$M %*% t(params$C) + 0.5 * (params$S*params$S) %*% t(params$C * params$C))
# pred <- exp(XB.hat + M.hat %*% t(C.hat) + 0.5 * (S.hat*S.hat) %*% t(C.hat * C.hat))
# 
# Y.logit <- ifelse(data$Y == 0, 0, 1)
# 
# #----------------------------  Plot 
# 
# #--------------- ELBO 
# par(mfrow=c(1, 1))
# plot(out$elboPath[out$elboPath > 1*out$elboPath[[1]]],type='l', main = "ELBO path",ylab='ELBO',xlab='iter')
# plot(D.hat, params$D) ; abline(0,1)
# 
# #----------------- Estimations
# 
# par(mfrow=c(1, 2))
# plot(XD, XD.init, col='red', main = "Estimation de la logistique") ; abline(0,1)
# points(XD, XD.hat) ; abline(0,1)
# 
# 
# 
# plot(XB, XB.init,col='red', main = "Estimation des régresseurs") ; abline(0,1)
# points(XB, XB.hat) ; abline(0,1)
# 
# 
# par(mfrow=c(1,3))
# boxplot(VectorToMatrix(XD, n, p) ~ Y.logit, main = "True")
# boxplot(VectorToMatrix(XD.init, n, p) ~ Y.logit, main = "Initialisation")
# boxplot(VectorToMatrix(XD.hat, n, p) ~ Y.logit, main = "Estimation")
# 
# 
# par(mfrow=c(2,2))
# 
# plot(log(1 + data$Y[data$Y !=0]), log(1 + pred[data$Y != 0]), main = "Prediction") ; abline(0,1)
# plot(log(1 + data$Y[data$Y !=0]), log(1 + pred.init[data$Y != 0]), main = "Prediction") ; abline(0,1)
# 
# plot( data$Y[data$Y !=0], pred[data$Y != 0], main = "Prediction") ; abline(0,1)
# plot( data$Y[data$Y !=0], pred.init[data$Y != 0], main = "Prediction") ; abline(0,1)



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






