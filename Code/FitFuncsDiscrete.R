library(nloptr)

FRouter <- function(n1, n2, Pval, c1, c2, h1, h2) {
  z <- matrix(NA, length(n1), length(Pval))
  for(i in 1:length(n1)) {
    for(j in 1:length(Pval)) {
      z[i,j] <- n1[i]*(1 - exp(-FR2species(N1=n1[i], N2=n2, P=Pval[j], c1=c1, c2=c2, h1=h1, h2=h2, type="2")))
       #N1*(1 - exp(-x1))
    }
  }
  return(z)
}


FR2speciesCT = function(N1, N2, P, c1, c2, h1, h2, type) {
  temp = switch(type, "1"=c1*N1, "2"=c1*N1/(1 + h1*c1*N1 + h2*c2*N2), "ratio"=c1*N1/(P + h1*c1*N1 + h2*c2*N2), "donor"=c1*N1/P)
  return(temp)
}

FR2species = function(N1, N2, P, c1, c2, h1, h2, type) {
  temp = switch(type, "1"=c1*P, "2"=c1*P/(1 + c1*h1*N1 + c2*h2*N2))
  return(temp)
}



optimFunc = function(parvec, objFunc, type, lb.vec, ub.vec) {

  if(!exists("N1SE")) {N1SE = NULL}
  if(!exists("N2SE")) {N2SE = NULL}
  if(!exists("PSE"))  {PSE  = NULL}
  if(!exists("DPSE")) {DPSE = NULL}
  par0=parvec

  #optim = try(crs2lm(x0=parvec, fn=objFunc, lower=parvec-10, upper=parvec+10, maxeval=1e5, N1=N1, N2=N2, Pred=P, DP=DP, type=type), silent=TRUE)
  #print(optim) NLOPT_GN_CRS2_LM NLOPT_GN_CRS2_LM NLOPT_GN_DIRECT_L_RAND
  #my.options <- list("algorithm"="NLOPT_GN_MLSL", "maxeval"=1e5, "local_opts"="NLOPT_LN_SBPLX")
  
  my.options <- list("algorithm"="NLOPT_GN_DIRECT_L", "ftol_rel"=1e-8, "maxeval"=1e5)
  optim.nloptr = try(nloptr(x0=parvec, eval_f=objFunc, lb=parvec-10, ub=parvec+10, N1=N1, N2=N2, Pred=P, DP=DP, type=type, opts=my.options), silent=TRUE)

  my.options <- list("algorithm"="NLOPT_LN_BOBYQA", "ftol_rel"=1e-8, "maxeval"=1e7)
  optim.nloptr = try(nloptr(x0=parvec, eval_f=objFunc, N1=N1, N2=N2, Pred=P, DP=DP, type=type, opts=my.options), silent=TRUE)

  
  my.options<-list("algorithm"="NLOPT_LN_SBPLX  ", "ftol_rel"=1e-8, "maxeval"=1e7)
  optim2.nloptr = try(nloptr(x0=optim.nloptr$solution, eval_f=objFunc, N1=N1, N2=N2, Pred=P, DP=DP, type=type, opts=my.options), silent=TRUE)


  if(class(optim2.nloptr) != 'try-error') {
        optim2 <- list()
        optim2$value <- optim2.nloptr$objective
        optim2$par <- optim2.nloptr$solution
        names(optim2$par) = names(parvec)
        optim2$convergence <- 0
        optim2$termination_conditions <- optim2.nloptr$termination_conditions
        optim2$status <- optim2.nloptr$status
        optim2$message <- optim2.nloptr$message
        optim2$iteractions <- optim2.nloptr$iterations
  } else {
    optim2 <- list()
    optim2$val <- NA
    optim2$par <- rep(NA, length(parvec))
    optim2$convergence <- -1
    return(optim2)
  }

  return(optim2)

}

logistic <- function(x) {
  1/(1 + exp(-x))
}


##fit dynamical models
mvnObj = function(parvec, Pred, N1, N2, DP, type) {

  P = Pred
  #type I response
  if(type == "1") {
    r1 = exp(parvec[1])
    s1 = exp(parvec[2])
    r2 = exp(parvec[3])
    s2 = exp(parvec[4])
    c1 = exp(parvec[5])
    h1 = 0
    e1 = logistic(parvec[6])
    c2 = exp(parvec[7])
    h2 = 0
    e2 = logistic(parvec[8])
    u = exp(parvec[9])
    sd1 = exp(parvec[10])
    sd2 = exp(parvec[11])
    sd3 = exp(parvec[12])
  }
  if(type == "2") {
    #type II response
    r1 = exp(parvec[1])
    s1 = exp(parvec[2])
    r2 = exp(parvec[3])
    s2 = exp(parvec[4])
    c1 = exp(parvec[5])
    h1 = exp(parvec[6])
    e1 = logistic(parvec[7])
    c2 = exp(parvec[8])
    h2 = exp(parvec[9])
    e2 = logistic(parvec[10])
    u = exp(parvec[11])
    sd1 = exp(parvec[12])
    sd2 = exp(parvec[13])
    sd3 = exp(parvec[14])
 }

  if(any(!is.finite(exp(parvec)))) {return(1e6)}

  Ppred  <- vector('numeric', length(P))
  N1pred <- vector('numeric', length(P))
  N2pred <- vector('numeric', length(P))


  x1 <- FR2species(N1=N1, N2=N2, P=P, c1=c1, c2=c2, h1=h1, h2=h2, type=type)
  x2 <- FR2species(N1=N2, N2=N1, P=P, c1=c2, c2=c1, h1=h2, h2=h1, type=type)

  N1prey   <- N1*(1 - exp(-x1))
  N2prey   <- N2*(1 - exp(-x2))
  
  Ppred   <- e1*N1prey + e2*N2prey + P*exp(-u)  
  N1pred  <- r1*N1*exp(-s1*N1 - x1)
  N2pred  <- r2*N2*exp(-s2*N2 - x2)

  if(any(is.nan(c(Ppred, N1pred, N2pred)))) {return(1e6)}
  abundLL = sum(dnorm(log(P[-1]), log(Ppred[-length(Ppred)]), sd1, log=TRUE) + dnorm(log(N1[-1]), log(N1pred[-length(N1pred)]), sd2, log=TRUE) + dnorm(log(N2[-1]), log(N2pred[-length(N2pred)]), sd3, log=TRUE))

  return(-abundLL)

}



##fit dynamical models coupled with diet data
mvnDietObj = function(parvec, Pred, N1, N2, DP, type) {

  P = Pred
  #type I response
  if(type == "1") {
    r1 <- exp(parvec[1])
    s1 <- exp(parvec[2])
    r2 <- exp(parvec[3])
    s2 <- exp(parvec[4])
    c1 <- exp(parvec[5])
    h1 <- 0
    e1 <- logistic(parvec[6])
    c2 <- exp(parvec[7])
    h2 <- 0
    e2 <- logistic(parvec[8])
    u  <- exp(parvec[9])
    sd1 <- exp(parvec[10])
    sd2 <- exp(parvec[11])
    sd3 <- exp(parvec[12])
    sdlogit <- exp(parvec[13])
  }
  if(type == "2") {
    #type II response
    r1 = exp(parvec[1])
    s1 = exp(parvec[2])
    r2 = exp(parvec[3])
    s2 = exp(parvec[4])
    c1 = exp(parvec[5])
    h1 = exp(parvec[6])
    e1 = logistic(parvec[7])
    c2 = exp(parvec[8])
    h2 = exp(parvec[9])
    e2 = logistic(parvec[10])
    u = exp(parvec[11])
    sd1 <- exp(parvec[12])
    sd2 <- exp(parvec[13])
    sd3 <- exp(parvec[14])
    sdlogit = exp(parvec[15])
  }

  if(any(!is.finite(exp(parvec)))) {return(1e6)}
  
  x1        <- FR2species(N1=N1, N2=N2, P=P, c1=c1, c2=c2, h1=h1, h2=h2, type=type)
  x2        <- FR2species(N1=N2, N2=N1, P=P, c1=c2, c2=c1, h1=h2, h2=h1, type=type)

  N1prey    <- N1*(1 - exp(-x1))
  N2prey    <- N2*(1 - exp(-x2))
  mu        <- (N1prey)/(N1prey + N2prey)

  Ppred   <- e1*N1prey + e2*N2prey + P*exp(-u)
  N1pred  <- r1*N1*exp(-s1*N1-x1)
  N2pred  <- r2*N2*exp(-s2*N2-x2)

  if(any(is.nan(c(Ppred, N1pred, N2pred, x1, x2)))) {return(1e6)}

  logitLL   <- dnorm(x=log(DP/(1-DP)), mean=log((mu)/(1-mu)), sd=sdlogit, log=TRUE)
  
  abundLL <- sum(dnorm(log(P[-1]), log(Ppred[-length(Ppred)]), sd1, log=TRUE) + dnorm(log(N1[-1]), log(N1pred[-length(N1pred)]), sd2, log=TRUE) + dnorm(log(N2[-1]), log(N2pred[-length(N2pred)]), sd3, log=TRUE))
  
  return( -sum(abundLL) - sum(logitLL) )

}

