#Simulate 1 predator 2 prey dynamics


sim.func = function(years, truePar, type, N0=NULL) {

  attach(truePar, warn.conflicts=F)
  P       <- vector('numeric', length=years)
  N1prey  <- vector('numeric', length=years)
  N2prey  <- vector('numeric', length=years)
  N1 <- vector('numeric', length=years)
  N2 <- vector('numeric', length=years)

  if(is.null(N0)) {
  if(type != "ratio" | type != "donor") {
    P[1]  <- 380
    N1[1] <- 500
    N2[1] <- 550
  } else {
    P[1]  <- 140
    N1[1] <- 265
    N2[1] <- 244
  }
  } else { 
    P[1]  <- N0[1]
    N1[1] <- N0[2]
    N2[1] <- N0[3]
  }

  for(i in 2:length(P)) {

    x1 <- FR2species(N1=N1[i-1], N2=N2[i-1], P=P[i-1], c1=c1[[1]], c2=c2[[1]], h1=h1[[1]], h2=h2[[1]], type=type)[[1]]
    x2 <- FR2species(N1=N2[i-1], N2=N1[i-1], P=P[i-1], c1=c2[[1]], c2=c1[[1]], h1=h2[[1]], h2=h1[[1]], type=type)[[1]]
    
    N1prey[i-1] <- N1[i-1]*(1 - exp(-x1))
    N2prey[i-1] <- N2[i-1]*(1 - exp(-x2))

    P[i]  <- (e1[[1]]*N1prey[i-1] + e2[[1]]*N2prey[i-1] + P[i-1]*exp(-u[[1]]))*exp(rnorm(1, 0, sigma[[1]]))
    N1[i] <- r1[[1]]*N1[i-1]*exp(-s1[[1]]*N1[i-1]-x1)*exp(rnorm(1, 0, sigma[[1]]))
    N2[i] <- r2[[1]]*N2[i-1]*exp(-s2[[1]]*N2[i-1]-x2)*exp(rnorm(1, 0, sigma[[1]]))

  }

  x1 <- FR2species(N1=N1[i], N2=N2[i], P=P[i], c1=c1[[1]], c2=c2[[1]], h1=h1[[1]], h2=h2[[1]], type=type)[[1]]
  x2 <- FR2species(N1=N2[i], N2=N1[i], P=P[i], c1=c2[[1]], c2=c1[[1]], h1=h2[[1]], h2=h1[[1]], type=type)[[1]]
    
  N1prey[i] <- N1[i]*(1 - exp(-x1))
  N2prey[i] <- N2[i]*(1 - exp(-x2))
     
  detach(truePar)

  return(list(P, N1, N2, N1prey, N2prey))

}


genData = function(samplesize, truePar, plotit=TRUE, type="2") {

  ##parameter definitions, from carpenter et al. where possible
  sim.out = sim.func(samplesize, truePar, type)

  P   <- sim.out[[1]]
  N1  <- sim.out[[2]]
  N2  <- sim.out[[3]]
  N1prey <- sim.out[[4]]
  N2prey <- sim.out[[5]]

  p1 <- N1prey/(N1prey+N2prey)
  p2 <- N2prey/(N1prey+N2prey)

  DP <- p1
  
  if(plotit==TRUE) {
    par(mfrow=c(1,2))
    plot(P, type='l', ylim=range(c(P,N1,N2), na.rm=T), lwd=2, xlab="Years", ylab="Abundances")
    lines(N1, col="cornflowerblue", lwd=2)
    lines(N2, col="red", lwd=2)

    plot(p1, type='l', lwd=2, col="cornflowerblue", xlab="Years", ylab="Diet proportion", ylim=c(0,1))
    lines(p2, lwd=2, col="red")
  }

  return(list(P=P, N1=N1, N2=N2, DP=DP))

}

