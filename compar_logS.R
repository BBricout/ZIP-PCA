rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

setwd("~/Documents/ZIP-PCA")

library(Rcpp)
library(PLNmodels)
library(missForest)

#--- Codes BB

source("codesBB/FunctionsBB.R")
source("codesBB/UtilsBB.R")


#--- Données

load(file='ourfile.Rdata')

n <- nrow(data$Y)
p <- ncol(data$Y)
q <- ncol(true$eStep$M)
d <- ncol(data$X)

tolXi <- 1e-04 ; tolS <- 1e-04

lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))

config <- PLNPCA_param()$config_optim
# config$algorithm <- "MMA" # Par défaut dans config c'est "CCSAQ"
config$lower_bounds <- lb

out.ref <- Miss.ZIPPCA(Y = data$Y.na, X = data$X, q, config = config)








