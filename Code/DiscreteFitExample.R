source('SimDynDiscrete.R')
source('FitFuncsDiscrete.R')
graphics.off()

#calculate SE's
batchFlag = !is.na(Sys.getenv("R_BATCH", NA))

plotit      <- FALSE
simlength   <- 5e2

args=(commandArgs(TRUE))

if(length(args)==0) {
    ##supply default values if not running in batch mode
    plotit  <- FALSE
    iters   <- 1e2
    samplesize <- 10
    type    <- '2'
    scale   <- 2
} else {
    for(i in 1:length(args)) {
      eval(parse(text=args[[i]]))
      plotit <- F
    }
}
cat(samplesize, type, '\n')

t1  <- data.frame(r1=1.8, s1=0.001, r2=1.8, s2=0.001, c1=0.0002, h1=0, e1=0.6, c2=0.00024, h2=0, e2=0.6, u=0.1, sigma=0.05*scale)
t2  <- data.frame(r1=1.8, s1=0.001, r2=1.8, s2=0.001, c1=0.001,  h1=5.0, e1=0.6, c2=0.003, h2=5.0, e2=0.6, u=0.1, sigma=0.05*scale)

truePar <- switch(type, "1"=t1, "2"=t2)

m1LL <- matrix(NA, iters, 2)
m2LL <- matrix(NA, iters, 2)
m1BIC <- matrix(NA, iters, 2)
m2BIC <- matrix(NA, iters, 2)

m1est   <- matrix(NA, iters, 14)
m2est   <- matrix(NA, iters, 15)

for(n in 1:iters) {

    set.seed(n)

    cat("iteration", n, '\n')

    out <- genData(samplesize=simlength, truePar=truePar, plotit=plotit, type=type)
    while(any(c(out$P, out$N1, out$N2) < 1)) {
         out <- genData(samplesize=simlength, truePar=truePar, plotit=plotit, type=type)
    }

    ind <- (simlength-samplesize+1):simlength
    P   <- out$P[ind]  
    N1  <- out$N1[ind] 
    N2  <- out$N2[ind] 
    DP  <- out$DP[ind] 
    p1  <- DP + rnorm(length(DP), 0, 0.001)
    if(any(p1 < 0)) {p1[which(p1 < 0)] <- 0}
    if(any(p1 > 1)) {p1[which(p1 > 1)] <- 1}
    DP  <- p1
    

    if(plotit) {
       # X11()
        n1vec   <- seq(min(N1), max(N1), length.out=50)
        n2vec   <- seq(min(N2), max(N2), length.out=50)
        pvec    <- seq(min(P), max(P), length.out=50)
        z       <- FRouter(n1vec, mean(N2), Pval=pvec, c1=t2$c1, c2=t2$c1, h1=t2$h1, h2=t2$h2)

       # persp(x=n1vec, y=pvec, z=z, lwd=2, theta=-0)
    }

    

    CV  <-  sd(P)/mean(P)
    if(!batchFlag) { print(sd(P)/mean(P)) 
        print(cor(N1, N2))
    }

    if(T) {
    ##Fit time series models, abundances
        #these are the starting values
    if(type =='1') {
        par1 <- c(r1=log(t1$r1), s1=log(t1$s1), r2=log(t1$r2), s2=log(t1$s2), c1=log(t1$c1), e1=log(t1$e1), c2=log(t1$c2), e2=log(t1$e2), u=log(t1$u), sd1=log(t1$sigma), sd2=log(0.5*t1$sigma), sd3=log(0.5*t1$sigma))
        par2 <- c(r1=log(t1$r1), s1=log(t1$s1), r2=log(t1$r2), s2=log(t1$s2), c1=log(t1$c1), h1=-10, e1=log(t1$e1), c2=log(t1$c2), h2=-10, e2=log(t1$e2), u=log(t1$u), sd1=log(t1$sigma), sd2=log(0.5*t1$sigma), sd3=log(0.5*t1$sigma))
    } else {
        par1 <- c(r1=log(t2$r1), s1=log(t2$s1), r2=log(t2$r2), s2=log(t2$s2), c1=log(t2$c1), e1=log(t2$e1), c2=log(t2$c2), e2=log(t2$e2), u=log(t2$u), sd1=log(t2$sigma), sd2=log(0.5*t2$sigma), sd3=log(0.5*t2$sigma))
        par2 <- c(r1=log(t2$r1), s1=log(t2$s1), r2=log(t2$r2), s2=log(t2$s2), c1=log(t2$c1), h1=log(t2$h1), e1=log(t2$e1), c2=log(t2$c2), h2=log(t2$h2), e2=log(t2$e2), u=log(t2$u), sd1=log(t2$sigma), sd2=log(0.5*t2$sigma), sd3=log(0.5*t2$sigma))
    }
        #now fit the models

    lb.vec  <- c(r1=-Inf, s1=-Inf, r2=-Inf, s2=-Inf, c1=-Inf, e1=-10, c2=-Inf, e2=-Inf, u=-Inf, sd1=-Inf, sd2=-Inf, sd3=-Inf)
    ub.vec  <- c(r1=Inf, s1=Inf, r2=Inf, s2=Inf, c1=Inf, e1=0, c2=Inf, e2=0, u=Inf, sd1=Inf, sd2=Inf, sd3=Inf)
    mvnType1GSA     <- optimFunc(parvec=par1, objFunc=mvnObj, type="1", lb.vec=lb.vec, ub.vec=ub.vec)
    lb.vec    <- c(r1=-Inf, s1=-Inf, r2=-Inf, s2=-Inf, c1=-Inf, h1=-Inf, e1=-Inf, c2=-Inf, h2=-Inf, e2=-Inf, u=-Inf, sd1=-Inf, sd2=-Inf, sd3=-Inf)
    ub.vec    <- c(r1=Inf, s1=Inf, r2=Inf, s2=Inf, c1=Inf, h1=Inf, e1=0, c2=Inf, h2=Inf, e2=0, u=Inf, sd1=Inf, sd2=Inf, sd3=Inf)
    mvnType2GSA     <- optimFunc(parvec=par2, objFunc=mvnObj, type="2", lb.vec=lb.vec, ub.vec=ub.vec)
   
        #check convergence
    #if(mvnType1GSA$convergence != 0 | mvnType2GSA$convergence != 0 || mvnType1GSA$value < mvnType2GSA$value) {   
    if(mvnType1GSA$convergence != 0 | mvnType2GSA$convergence != 0) {   
        cat('convergence failed', mvnType1GSA$value, mvnType2GSA$value, '\n')
        mvnType1GSA$value <- NA
        mvnType2GSA$value <- NA
    }
    LLvec       <- c(mvnType1GSA$value, mvnType2GSA$value) 
    n.sample    <- 3*length(p1[-1])
    k1          <- length(mvnType1GSA$par)
    k2          <- length(mvnType2GSA$par)
    kvec        <- c(k1,k2) 
    BICvec      <- 2*LLvec + kvec*log(n.sample)
  
    if(!batchFlag) { 
      cat("\nLL :", LLvec, '\n')
      cat("BIC:", BICvec - min(BICvec, na.rm=T), '\n')
    }

    m1LL[n,]    <- LLvec
    m1BIC[n,]   <- BICvec
    m1est[n,]   <- mvnType2GSA$par
    }   

    
    ##Fit joint diet abundance time series models
        #set starting values
    if(type =='1') {
        par1    <- c(par1[1:12], sdlogit=0) 
        par2    <- c(par2[1:5], h1=-10, par2[7:8], h2=-10, par2[10:14], sdlogit=0)
    } else {
        par1    <- c(par2[1:5], par2[7:8], par2[10:14], sdlogit=0) 
        par2    <- c(par2, sdlogit=0)
    }
        
        #fit models
   # lb.vec  <- c(r1=-1, s1=-100, r2=-1, s2=-100, c1=-20, e1=-10, c2=-20, e2=-10, u=-10, sd1=-10, sd2=-10, sd3=-10, sdlogit=-Inf)
    #ub.ve c <- c(r1=5, s1=0, r2=5, s2=0, c1=0, e1=0, c2=0, e2=0, u=0, sd1=10, sd2=10, sd3=10, sdlogit=Inf)
    lb.vec  <- c(r1=-Inf, s1=-Inf, r2=-Inf, s2=-Inf, c1=-Inf, e1=-10, c2=-Inf, e2=-Inf, u=-Inf, sd1=-Inf, sd2=-Inf, sd3=-Inf, sdlogit=-Inf)
    ub.vec  <- c(r1=Inf, s1=Inf, r2=Inf, s2=Inf, c1=Inf, e1=0, c2=Inf, e2=0, u=Inf, sd1=Inf, sd2=Inf, sd3=Inf, sdlogit=Inf)
    jointmvnType1GSA <- optimFunc(parvec=par1, objFunc=mvnDietObj, type="1", lb.vec=lb.vec, ub.vec=ub.vec)
    #lb.vec    <- c(r1=-1, s1=-100, r2=-1, s2=-100, c1=-20, h1=-10, e1=-10, c2=-20, h2=-10, e2=-10, u=-10, sd1=-10, sd2=-10, sd3=-10, sdlogit=-Inf)
    #ub.vec    <- c(r1=5, s1=0, r2=5, s2=0, c1=0, h1=5, e1=0, c2=0, h2=5, e2=0, u=0, sd1=10, sd2=10, sd3=10, sdlogit=Inf)
    lb.vec    <- c(r1=-Inf, s1=-Inf, r2=-Inf, s2=-Inf, c1=-Inf, h1=-Inf, e1=-Inf, c2=-Inf, h2=-Inf, e2=-Inf, u=-Inf, sd1=-Inf, sd2=-Inf, sd3=-Inf, sdlogit=-Inf)
    ub.vec    <- c(r1=Inf, s1=Inf, r2=Inf, s2=Inf, c1=Inf, h1=Inf, e1=0, c2=Inf, h2=Inf, e2=0, u=Inf, sd1=Inf, sd2=Inf, sd3=Inf, sdlogit=Inf)
    jointmvnType2GSA <- optimFunc(parvec=par2, objFunc=mvnDietObj, type="2", lb.vec=lb.vec, ub.vec=ub.vec)
    
         #check convergence
    #if(jointmvnType1GSA$convergence != 0 | jointmvnType2GSA$convergence != 0 || jointmvnType1GSA$value < jointmvnType2GSA$value) {   
    if(jointmvnType1GSA$convergence != 0 | jointmvnType2GSA$convergence != 0) {   
        cat('convergence failed', jointmvnType1GSA$value, jointmvnType2GSA$value, '\n')
        jointmvnType1GSA$value <- NA
        jointmvnType2GSA$value <- NA
    }

    LLvec <- c(jointmvnType1GSA$value, jointmvnType2GSA$value)

    n.sample    <- 3*length(p1[-1]) + length(p1)
    k1          <- length(jointmvnType1GSA$par)
    k2          <- length(jointmvnType2GSA$par)
    kvec        <- c(k1, k2) 
    BICvec      <- 2*LLvec + kvec*log(n.sample)

    m2LL[n,]    <- LLvec
    m2BIC[n,]   <- BICvec
    m2est[n,]   <- jointmvnType2GSA$par
    
    if(!batchFlag) {
        cat("\nLL :", LLvec, '\n')
        cat("BIC:", BICvec - min(BICvec, na.rm=T), '\n')
    }
    


    graphics.off()

    if(batchFlag) {
        save.image(file=paste("SimOut/IntenseComm_Nobs", samplesize, "scale", scale, "_type",type, ".Rdata",sep=''))
    }
}

if(batchFlag) {
    q('no')
}