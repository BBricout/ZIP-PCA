# Test with/without miss

rm(list=ls()); par(mfrow=c(1, 1), pch=20); palette('R3')
# seed <- .Random.seed
library(nloptr)
source('../Functions/FunctionsZIP.R'); 
source('../Functions/FunctionsUtils.R'); 
source('../Functions/FunctionsZIPLNmiss.R')
source('../Functions/FunctionsZIPLNmissVec.R')
simDir <- '../../simulSR/'
resDir <- './'

# Data
n <- 100; d <- 5; p <- 10; q <- 2; seed <- 5; obs <- 1
simParmsFull <- paste0('-n', n, '-d', d, '-p', p, '-q', q, '-seed', seed)
simParms <- paste0(simParmsFull, '-obs', 100*obs)
simName <- paste0('ZiPLNsim', simParms)
simFile <- paste0(simDir, simName, '.Rdata')
load(simFile)

# Algo parms
iMiss <- 22; jMiss <- 3; ijMiss <- c(iMiss, jMiss)
tolXi <- 0; tolS <- 1e-6

# Data sets
dataFull <- dataMiss1 <- dataMissRow <- dataAbsRow <- data; 
dataMiss1$Omega[iMiss, jMiss] <- 0
dataMissRow$Omega[iMiss, ] <- 0
dataAbsRow$Y <- data$Y[-iMiss, ]
dataAbsRow$logFactY <- data$logFactY[-iMiss, ]
dataAbsRow$Omega <- data$Omega[-iMiss, ]
dataAbsRow$X <- data$X[-which(data$ij[, 1] %in% iMiss), ]
dataAbsRow$ij <- cbind(rep(1:(n-1), p), rep(1:p, each=(n-1)))
dataiMiss <- list(Yi=dataMissRow$Y[iMiss, ], Xi=dataMissRow$X[which(data$ij[, 1]==iMiss), ], 
                  Omegai=dataMissRow$Omega[iMiss, ], logFactYi=dataMissRow$logFactY[iMiss, ])

# Init
init <- InitZiPLN(data=dataMissRow, q=q, tolXi=tolXi)
mStep <- init$mStep
eStep <- VEstep(data=dataMissRow, mStep=mStep, eStep=init$eStep, tolXi=tolXi, tolS=tolS)
initAbsRowEstep <- init$eStep; initAbsRowEstep$M <- initAbsRowEstep$M[-iMiss, ]; 
initAbsRowEstep$S <- initAbsRowEstep$S[-iMiss, ]; initAbsRowEstep$xi <- initAbsRowEstep$xi[-iMiss, ]
eStepAbsRow <- VEstep(data=dataAbsRow, mStep=mStep, eStep=initAbsRowEstep, tolXi=tolXi, tolS=tolS)
eStepiMiss <- list(xii=eStep$xi[iMiss, ], 
                   mi=as.vector(eStep$M[iMiss, , drop=FALSE]), 
                   Si=as.vector(eStep$S[iMiss, , drop=FALSE]))

# ELBO
elboMissRow <- ELBO(data=dataMissRow, mStep=mStep, eStep=eStep) 
elboAbsRow <- ELBO(data=dataAbsRow, mStep=mStep, eStep=eStepAbsRow) 
elboImiss <- ELBOi(datai=dataiMiss, mStep=mStep, eStepi=eStepiMiss)
c(elboMissRow, elboAbsRow, elboImiss, elboMissRow - elboAbsRow - elboImiss)

# Compare
rbind(eStep$M[iMiss, ], eStep$S[iMiss, ])
rbind(eStep$xi[iMiss, ], 
      as.vector(plogis(data$X[which(data$ij[, 1] %in% iMiss), ]%*%mStep$gamma)))
