library(runjags)

source('WhaleMCMC.jags')

Diets       <- read.csv(file="../Data/WhaleDietProportions.csv", header=T)
#Abundances  <- read.csv(file="../Data/WhaleAbundances.csv", header=T)
#WhaleSE     <- read.csv(file="../Data/WhaleAbundancesSE.csv", header=T)[,3]^2/Abundances[,3]^2
Abundances  <- read.csv(file="../Data/WhaleAbundancesSE.csv", header=T)
WhaleSE     <- read.csv(file="../Data/WhaleAbundancesSE.csv", header=T)[,8]#^2/Abundances[,3]^2
FishSE      <- read.csv(file="../Data/WhaleAbundancesSE.csv", header=T)[,9]
ZooSE       <- read.csv(file="../Data/WhaleAbundancesSE.csv", header=T)[,10]


bm1 = 137*1*10^3 #biomass of a fish adjusted for digestability
bm2 = 15.5*0.93*10^-3
tstep <- c(NA, 1, 2, 5, 1, 1)


##Fit time series models
P     <<- Abundances[,3]
Pred  <<- P
N1    <<- Abundances[,5]
N2    <<- Abundances[,6]
DP    <<- Diets[,9]/100
N	    <<- Diets$N

dt    <- 0.01
tseq  <- seq(0, 1, by=dt)
na.vec <- which(is.na(N1))
datTable <- list(N1obs=N1[-na.vec], N2obs=N2[-na.vec], P=P[-na.vec], PSE=WhaleSE[-na.vec], N1SE=FishSE[-na.vec], N2SE=ZooSE[-na.vec], DP=DP[-na.vec], DP.sd=Diets[-na.vec,10]/100, DPlogit= log(DP[-na.vec]/(1-DP[-na.vec])), N=6, tstep=tstep, tseq=tseq, dt=dt, lambda=log(2)/(14/365), bm1=bm1, bm2=bm2, DP.samplesize=length(DP[-na.vec]))

par.est  <- c('r1', 'r2', 'c1', 'c2', 'procE1', 'procE2', 'tauDraw', 'deviance') 

abund.ss  <- run.jags(ssModel, burnin=5e4, sample=1e5, n.chains=4, thin=10, adapt=1e3, monitor=par.est, data=datTable, method="parallel")

int.ss  <- run.jags(ssIntModel, burnin=5e4, sample=5e4, n.chains=4, thin=10, adapt=1e3, monitor=par.est, data=datTable, method="parallel")

print(int.ss)

#save.image(file="WhaleComparisonTypeIbm1.Rdata")
save.image(file="WhaleComparisonTypeI_slow.Rdata")

