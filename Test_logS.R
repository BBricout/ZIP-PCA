rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')

setwd("~/Documents/ZIP-PCA") # A changer
library(Rcpp)
library(PLNmodels)

source("codesBB/UtilsBB.R")
source("Data_create.R")
sourceCpp("src/optim_rank_ZIP_log.cpp")


init <- InitZiPLN(data,q) 
params <- list(B = as.matrix(init$mStep$beta), 
               D = as.matrix(init$mStep$gamma), 
               C = as.matrix(init$mStep$C), 
               M = as.matrix(init$eStep$M), 
               S = as.matrix(init$eStep$S))

paramsLogS <- list(B = as.matrix(init$mStep$beta), 
               D = as.matrix(init$mStep$gamma), 
               C = as.matrix(init$mStep$C), 
               M = as.matrix(init$eStep$M), 
               logS = log(as.matrix(init$eStep$S)))
config <- PLNPCA_param()$config_optim ; tolXi <- 1e-04


ElboB(data, params, tolXi)$objective
ElboBLogS(data, paramsLogS, tolXi)$objective


