
rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
library(Rcpp)
library(PLNmodels)

# seed <- .Random.seed
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIPLN.R')

source('codesSR/Functions/FunctionsZIP.R')

simDir <- 'SimulationsBB/datasim/'
# Parms: many small sims
n <- 100; d <- 5; p <- 10; q <- 2
baseSimName <- 'ZiPLNsim'
seedList <- 1:10; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)

# # Parms: one big sim
# n <- 500; d <- 20; p <- 30; q <- 5
# seedList <- 1; seedNb <- length(seedList)
# obsList <- c(0.6); obsNb <- length(obsList)

#X0 <- NULL
# Same X for all
seedX <- 1; set.seed(seedX)
X0 <- matrix(rnorm(n*p*d), n*p, d); X0[, 1] <- 1; 
baseSimName <- paste0(baseSimName, '-sameX', seedX)




#############################################################
set.seed(1)
seed = 1
simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
simNameFull <- paste0(baseSimName, simParmsFull)
simFileFull <- paste0(simDir, simNameFull, '-noMiss.Rdata')
sim <- SimZiPLN(n=n, p=p, d=d, q=q, X=X0)
data <- list(X=sim$X, Y=sim$Y,  ij=sim$ij, logFactY=lgamma(sim$Y+1))
true <- list(mStep=list(gamma=sim$gamma, beta=sim$beta, C=sim$C),
             eStep=list(M=sim$W, S=matrix(1e-4, n, q), xi=matrix(plogis(sim$X%*%sim$gamma), n, p)),
             latent=list(U=sim$U, W=sim$W, Z=sim$Z, Yall=sim$Yall))

data$Omega <- matrix(1, n, p)
data$R <- data$Omega


simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
simNameFull <- paste0(baseSimName, simParmsFull)
simFileFull <- paste0(simDir, simNameFull, '-noMiss.Rdata')
load(file=simFileFull)

Y <- data$Y
X <- data$X

source("FunctionsBB.R")

# config <- list(
#   algorithm = "CCSAQ",
#   backend = "nlopt",
#   maxeval = 10000,
#   ftol_rel = 1e-08,
#   xtol_rel = 1e-06,
#   ftol_abs = 0,
#   xtol_abs = 0,
#   maxtime = -1,
#   trace = 1,
#   lower_bound_S = matrix(1e-6, nrow = n, ncol = q) # Exemple de borne inférieure
# )
# config$maxeval <- 3

params <- Init_ZIP(Y, X, q)
params$S <- exp(params$S)

# start_timeB <- Sys.time()
outB <- Miss.ZIPPCA(data$Y, data$X, q, tolXi = 1e-4)
# end_timeB <- Sys.time()
# execution_time_realB <- end_timeB - start_timeB

outB$elboPath

plot(log(1 + data$Y[data$Y != 0]), log(1 + outB$pred$A[data$Y != 0])) ; abline(0,1)

plot(outB$params.init$B, outB$mStep$beta) ; abline(0,1)
plot(outB$params.init$D, outB$mStep$gamma) ; abline(0,1)
plot(outB$params.init$C, outB$mStep$C) ; abline(0,1)
plot(outB$params.init$M, outB$eStep$M) ; abline(0,1)
plot(outB$params.init$S, outB$eStep$S) ; abline(0,1)
