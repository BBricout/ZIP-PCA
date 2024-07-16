rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

setwd("~/Documents/ZIP-PCA")

library(Rcpp)
library(PLNmodels)
library(missForest)



#--- Codes SR
source('codesSR/Functions/FunctionsUtils.R')
source('codesSR/Functions/FunctionsZIP.R')
source("codesSR/Functions/FunctionsZIPLNmiss.R")
source("codesSR/Functions/FunctionsZIPLNmissVec.R")

#--- Codes BB

# sourceCpp("src/optim_rank_ZIP_new.cpp")
source("codesBB/FunctionsBB.R")
source("codesBB/UtilsBB.R")

#--- Données

load(file='ourfile.Rdata')

n <- nrow(data$Y)
p <- ncol(data$Y)
q <- ncol(true$eStep$M)
d <- ncol(data$X)

tolXi <- 1e-04 ; tolS <- 1e-06
lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))
config <- PLNPCA_param()$config_optim
config$lower_bounds <- lb


#--- Application aux jeux de données 


out.ref <- Miss.ZIPPCA(Y = data$Y, X = data$X, q=q, config = config)
out.log <- Miss.ZIPPCA.logS(Y=data$Y, X=data$X, q=q)
plot(out.ref$elboPath, col = "red")
points(out.log$elboPath)
out.ref$elbo
out.log$elbo



Y.test <- prodNA(data$Y, 0.1)

names(Init_ZIP(Y.test, data$X, q))
out.ref.na <- Miss.ZIPPCA(Y = Y.test, X = data$X, q=q, config=config)
out.log.na <- Miss.ZIPPCA.logS(Y=Y.test, X=data$X, q=q)




