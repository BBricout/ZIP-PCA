# Test with/without miss

rm(list=ls())
library(nloptr)
setwd('/home/robin/RECHERCHE/ECOLOGIE/B-Bricout/ZIP-PCA/codesSR/tests/')
source('../Functions/FunctionsZIP.R'); 
source('../Functions/FunctionsUtils.R'); 
source('../Functions/FunctionsZIPLNmiss.R')
source('../Functions/FunctionsZIPLNmissVec.R')
source('../Functions/FunctionsBB2SR.R')
setwd('/home/robin/RECHERCHE/ECOLOGIE/B-Bricout/ZIP-PCA/')
source('codesBB/FunctionsBB.R')
source('codesBB/UtilsBB.R')
setwd('/home/robin/RECHERCHE/ECOLOGIE/B-Bricout/ZIP-PCA/codesSR/tests/')
simDir <- '../../simulSR/'
resDir <- './'

# Data
n <- 100; p <- 10; d <- 5; q <- 2; seed <- 1; obs <- 0.5
load(paste0(simDir, 'ZiPLNsim-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed, '-obs', round(100*obs), '.Rdata'))

# Algo parms
iterMax <- 1e3; tolS <- 1e-6; tolXi <- 1e-4
lb <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))
lbAbs <- c(rep(-Inf, 2*d + q*(p+n-1)), rep(tolS, (n-1) * q))
algoName <- paste0('-tolS', tolS, '-tolXi', tolXi)
optsMMA <- list("algorithm"="NLOPT_LD_MMA", "xtol_rel"=1.0e-8)
optsLBGFS <- list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8)
optsCCSAQ <- list("algorithm"="NLOPT_LD_CCSAQ", "xtol_rel"=1.0e-8)

# Functions
Init2logS <- function(init){init$eStep$S <- log(init$eStep$S); return(init)}

# Fit SR, S
init <- InitZiPLN(data=data, q=q, tolXi=tolXi)
thetaPsi <- c(Mstep2Theta(init$mStep, n=n, d=d, p=p, q=q), Estep2Psi(init$eStep, n=n, d=d, p=p, q=q))
lb.sr <- c(rep(-Inf, 2*d + q*(p+n)), rep(tolS, n * q))
mma.sr <- try(nloptr(x0=thetaPsi,
                    eval_f=NegElboVecThetaPsi, eval_grad_f=NegElboGradVecThetaPsi, data=data, q=q, tolXi=tolXi, lb=lb, opts=optsLBGFS))
lbfgs.sr <- try(nloptr(x0=thetaPsi,
                    eval_f=NegElboVecThetaPsi, eval_grad_f=NegElboGradVecThetaPsi, data=data, q=q, tolXi=tolXi, lb=lb, opts=optsMMA))
ccsaq.sr <- nloptr(x0=thetaPsi,
                   eval_f=NegElboVecThetaPsi, eval_grad_f=NegElboGradVecThetaPsi, data=data, q=q, tolXi=tolXi, opts=optsCCSAQ)
gd.sr <- try(optim(par=thetaPsi,
                fn=ElboVecThetaPsi, gr=ElboGradVecThetaPsi, data=data, tolXi=tolXi,
                      q=q, method='L-BFGS-B', control=list(fnscale=-1, maxit=iterMax), lower=lb))
# Fit SR, logS
initLogS <-Init2logS(init)
thetaPsiLogS <- c(Mstep2Theta(initLogS$mStep, n=n, d=d, p=p, q=q), Estep2Psi(initLogS$eStep, n=n, d=d, p=p, q=q))
mmaLogS.sr <- nloptr(x0=thetaPsiLogS,
                  eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=data, q=q, tolXi=tolXi, opts=optsMMA)
lbfgsLogS.sr <- nloptr(x0=thetaPsiLogS,
                    eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=data, q=q, tolXi=tolXi, opts=optsLBGFS)
ccsaqLogS.sr <- nloptr(x0=thetaPsiLogS,
                    eval_f=NegElboLogSVecThetaPsi, eval_grad_f=NegElboLogSGradVecThetaPsi, data=data, q=q, tolXi=tolXi, opts=optsCCSAQ)
gdLogS.sr <- optim(par=thetaPsiLogS,
                fn=ElboLogSVecThetaPsi, gr=ElboLogSGradVecThetaPsi, data=data, tolXi=tolXi,
                q=q, method='BFGS', control=list(fnscale=-1, maxit=iterMax))

# Fit BB, S
data$R <- data$Omega; 
data$Y.na <- data$Y; data$Y.na[which(data$Omega==0)] <- NA
mma.bb <- try(nloptr(x0=thetaPsi,
                     eval_f=Bobj.neg, eval_grad_f=Bgrad.neg, data=data, tolXi=tolXi, lb=lb, opts=optsLBGFS))
lbfgs.bb <- try(nloptr(x0=thetaPsi,
                       eval_f=Bobj.neg, eval_grad_f=Bgrad.neg, data=data, tolXi=tolXi, lb=lb, opts=optsMMA))
ccsaq.bb <- nloptr(x0=thetaPsi,
                       eval_f=Bobj.neg, eval_grad_f=Bgrad.neg, data=data, tolXi=tolXi, opts=optsCCSAQ)
gd.bb <- try(optim(par=thetaPsi,
                   fn=Bobj, gr=Bgrad, data=data, tolXi=tolXi,
                   method='L-BFGS-B', control=list(fnscale=-1, maxit=iterMax), lower=lb))
# Fit BB, logS
mmaLogS.bb <- nloptr(x0=thetaPsiLogS,
                     eval_f=BobjLogS.neg, eval_grad_f=BgradLogS.neg, data=data, tolXi=tolXi, opts=optsMMA)
lbfgsLogS.bb <- nloptr(x0=thetaPsiLogS,
                       eval_f=BobjLogS.neg, eval_grad_f=BgradLogS.neg, data=data, tolXi=tolXi, opts=optsLBGFS)
ccsaqLogS.bb <- nloptr(x0=thetaPsiLogS,
                       eval_f=BobjLogS.neg, eval_grad_f=BgradLogS.neg, data=data, tolXi=tolXi, opts=optsCCSAQ)
gdLogS.bb <- optim(par=thetaPsiLogS,
                   fn=BobjLogS, gr=BgradLogS, data=data, tolXi=tolXi,
                   method='BFGS', control=list(fnscale=-1, maxit=iterMax))

res <- rbind(c(-mma.sr$objective, -lbfgs.sr$objective, -ccsaq.sr$objective, gd.sr$value), 
             c(-mmaLogS.sr$objective, -lbfgsLogS.sr$objective, -ccsaqLogS.sr$objective, gdLogS.sr$value), 
             c(-mma.bb$objective, -lbfgs.bb$objective, -ccsaq.bb$objective, gd.bb$value), 
             c(-mmaLogS.bb$objective, -lbfgsLogS.bb$objective, -ccsaqLogS.bb$objective, gdLogS.bb$value))
colnames(res) <- c('mma', 'lbfgs', 'ccsaq', 'optim')
rownames(res) <- c('SR', 'SRlogS', 'BB', 'BBlogS')
res
