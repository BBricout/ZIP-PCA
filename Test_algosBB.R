#-------------Tests des algos-------------------



## Préparation d'un jeu de données

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
library(Rcpp)
library(PLNmodels)

# seed <- .Random.seed
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIPLN.R')
source("FunctionsBB.R")
source('codesSR/Functions/FunctionsZIP.R')

simDir <- 'SimulationsBB/datasim/'
# Parms: many small sims
n <- 100; d <- 5; p <- 10; q <- 2
baseSimName <- 'ZiPLNsim'
seedList <- 1:10; seedNb <- length(seedList)
obsList <- c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5); obsNb <- length(obsList)

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
config <- PLNPCA_param()$config_optim ; tolXi <- 0
config$algorithm <- "MMA"
init <- InitZiPLN(data, q, tolXi = 0) 
params <- list(B = as.matrix(init$mStep$beta), 
               D = as.matrix(init$mStep$gamma), 
               C = as.matrix(init$mStep$C), 
               M = as.matrix(init$eStep$M), 
               S = as.matrix(init$eStep$S))

plot(init$mStep$beta, params$B) ; abline(0,1)
plot(init$mStep$gamma, params$D) ; abline(0,1)
plot(init$mStep$C, params$C) ; abline(0,1)
plot(init$eStep$M, params$M) ; abline(0,1)
plot(init$eStep$S, params$S) ; abline(0,1)

## optim_rank_ZIP

sourceCpp(file = "src/optim_rank_ZIP.cpp")

out_ref <- nlopt_optimize_ZIP(data, params, config, tolXi)
out_ref$objective_values[length(out_ref$objective_values)]
plot(out_ref$objective_values)

# vem <- VemZiPLN(data=data, init=init, iterMax=5e3)

range(out_ref$S)
XB <- data$X%*%true$mStep$beta
XB_ref <- data$X%*%out_ref$B
XB_vem <- data$X%*%vem$mStep$beta
plot(XB_ref, XB) ; abline(0,1)
plot(XB_vem, XB) ; abline(0,1)

cbind(ELBOSophie(data = data, init$mStep,  init$eStep), unlist(ElboB(data,params, tolXi)))
Selbo <- ELBO(data, init$mStep, init$eStep)


## optim avec lowerbound

lb1 <- rep(-Inf, 2*d + q*(p+n) + 2)
lb2 <- rep(1e-04, n * q)
lb <- c(lb1, lb2)
sourceCpp("src/optim_rank_ZIP_lb.cpp")

config$lower_bounds <- lb
out_lb <- nlopt_optimize_ZIP_lb(data, params, config, tolXi)
out_lb$objective_values[length(out_lb$objective_values)]
out_lb$objective_values











