# Sim and fit ZI-PLN

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
#save(data, true, file=simFileFull)

data$Omega <- ifelse(is.na(data$Y), 0, 1)
data$R <- data$Omega


#############################################################
################ INIT ajustement SR
init <- InitZiPLN(data) #! Fichiers différnets pour miss et pas miss
###################################################



params_init <- list(B = as.matrix(init$mStep$beta), 
               D = as.matrix(init$mStep$gamma), 
               C = as.matrix(init$mStep$C), 
               M = as.matrix(init$eStep$M), 
               S = as.matrix(init$eStep$S))

data <- list(Y = as.matrix(data$Y), 
              R = as.matrix(ifelse(is.na(data$Y), 0, 1)), 
              X = as.matrix(data$X),ij=sim$ij,  logFactY=lgamma(sim$Y+1))

config <- PLNPCA_param()$config_optim

############################################# 
########### COMPAR ELBO 

source("FunctionsBB.R")

outB <- Miss.ZIPPCA(data$Y, data$X, q, params = params_init)

params = c(outB$mStep,outB$eStep)
params$B <- params$beta
params$D <- params$gamma


eStep = outB$eStep
mStep = outB$mStep
ELBOSophie(data = data, mStep,  eStep )

cbind(ELBOSophie(data = data, mStep,  eStep), unlist(ElboB(data,params)))


mStep = outB$mStep 
eStep = outB$eStep
#ElboSfun <- ELBO(data = data, mStep = outB$mStep, eStep = outB$eStep)

ElboSfun <- ELBOSophie(data = data, mStep = outB$mStep, eStep = outB$eStep)
ElboSfun[6]













