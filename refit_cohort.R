############# refit algorithm by Xin Wang#######
# refit based on a given group information
# year and age should be ordered increasingly
# group is the predefined group information, the order is the same as index (year or age)
# model can be year or age. If year is specified, group is for year
# since equal size the weighted square is not necessary 



refit_cohort <- function(year, age, y, x, group, model = "year")
{
  
  n0 <- length(y) # number of total observation
  ncx <- ncol(x)
  
  uniq_year <- unique(year) 
  uniq_age <- unique(age)
  
  nyear <- length(uniq_year)
  nage <- length(uniq_age)
  ncoh <- nyear + nage - 1 # number of cohort effects
  
  if(model == "year")
  {
    nobs <- nyear ## number of individuals
    nrep <- nage # number of observations for each indiviual
    indexy <- year
    uniq_index <- uniq_year
  }
 
  if(model == "age")
  {
    nobs <- nage
    nrep <- nyear
    indexy <- age
    uniq_index <- uniq_age
  }
  
  cohort <- as.integer(factor(year - age))
  Ip <- diag(1,ncx,ncx)
  
  Zc <- matrix(0, n0, ncoh) # cohort matrix
  Zc[cbind(1:n0, cohort)] <- 1
  
  ## 
  Xm <- matrix(0, n0, nobs*ncx)
  for(i in 1:nobs)
  {
    indexi <- indexy == uniq_index[i]
    Xm[indexi, (ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniq_index[i],]
  }
  
  ng <- length(unique(group))
  W <- matrix(0, nobs, ng)
  W[cbind(1:nobs,group)] <- 1
  W <- W %x% diag(1,ncx)
  Ux <- cbind(Xm%*%W, Zc) ## new X with group information

  ncu <- ncol(Ux)
  
  if(ng < nobs)
  {
    H1 <- matrix(0, nrow = 1, ncol = ncu) # constraint including beta 
    H1[,(ng*ncx+1):ncu] <- 1
    
    U1 <- matrix(0, ncu + 1,  ncol = ncu + 1) # linear model with constraint
    U1[1:ncu, 1:ncu] <- t(Ux)%*%Ux
    
    U1[ncu + 1, 1:ncu] <- H1
    U1[1:ncu, ncu+1] <- t(H1)
    thetaest <- (solve(U1)%*%c(t(Ux)%*%y,0))[1:ncu]
    
  }
  
  if(ng == nobs)
  {
    # three constraints
    H1 <- matrix(0, nrow = 3, ncol = ncu) # constraint including beta 
    H1[1,(ng*ncx+1):ncu] <- 1
    H1[2,(ng*ncx+1):ncu] <- 1:ncoh
    H1[3,(ng*ncx+1):ncu] <- (1:ncoh)^2
    
    U1 <- matrix(0, ncu + 3,  ncol = ncu + 3) # linear model with constraint
    U1[1:ncu, 1:ncu] <- t(Ux)%*%Ux
    
    U1[ncu + 1:3, 1:ncu] <- H1
    U1[1:ncu, ncu+1:3] <- t(H1)
    thetaest <- (solve(U1)%*%c(t(Ux)%*%y,0,0,0))[1:ncu]
  }
  

  # U_inv <- solve(t(Ux)%*%Ux + t(H1)%*%H1)
  # thetaest <- (U_inv - U_inv%*%t(H1)%*% solve(H1%*%U_inv%*%t(H1))%*%H1%*%U_inv)%*%t(Ux)%*%y # solution 

  
  betaest <- matrix(thetaest[1:(ng*ncx)],ncol = ncx, byrow=TRUE)
  etaest <- thetaest[-(1:(ng*ncx))]

  estimates <- Ux %*% thetaest
  
  sig2 <- sum((y - estimates)^2)/n0
  
  curve_est <- Xm%*%W%*%thetaest[1:(ng*ncx)]

  out <- list(betaest = betaest, etaest = etaest, estimates = estimates, curve = curve_est, sig2 = sig2)
  return(out)
}


##### refit without cohort effect ####

refit_group <- function(year, age, y, x, group, model = "year")
{
  
  n0 <- length(y) # number of total observation
  ncx <- ncol(x)
  
  uniq_year <- unique(year) 
  uniq_age <- unique(age)
  
  nyear <- length(uniq_year)
  nage <- length(uniq_age)
  ncoh <- nyear + nage - 1 # number of cohort effects
  
  if(model == "year")
  {
    nobs <- nyear ## number of individuals
    nrep <- nage # number of observations for each indiviual
    indexy <- year
    uniq_index <- uniq_year
  }
  
  if(model == "age")
  {
    nobs <- nage
    nrep <- nyear
    indexy <- age
    uniq_index <- uniq_age
  }
  
 
  Ip <- diag(1,ncx,ncx)
  
  ## 
  Xm <- matrix(0, n0, nobs*ncx)
  for(i in 1:nobs)
  {
    indexi <- indexy == uniq_index[i]
    Xm[indexi, (ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniq_index[i],]
  }
  
  ng <- length(unique(group))
  W <- matrix(0, nobs, ng)
  W[cbind(1:nobs,group)] <- 1
  W <- W %x% diag(1,ncx)
  Ux <- Xm%*%W  ## new X with group information
  
  betaest <- solve(t(Ux)%*%Ux)%*%t(Ux)%*%y
  
  estimates <- Ux %*% betaest
  sig2 <- sum((y - estimates)^2)/n0


  out <- list(betaest = betaest, estimates = estimates, sig2 = sig2)
  return(out)
}


##### refit with given group for both individual and cohort ####
refit_cohort2 <- function(year, age, y, x, group.individual, group.cohort, model = "year")
{
  
  n0 <- length(y) # number of total observation
  ncx <- ncol(x)
  
  uniq_year <- unique(year) 
  uniq_age <- unique(age)
  
  nyear <- length(uniq_year)
  nage <- length(uniq_age)
  ncoh <- nyear + nage - 1 # number of cohort effects
  
  if(model == "year")
  {
    nobs <- nyear ## number of individuals
    nrep <- nage # number of observations for each indiviual
    indexy <- year
    uniq_index <- uniq_year
  }
  
  if(model == "age")
  {
    nobs <- nage
    nrep <- nyear
    indexy <- age
    uniq_index <- uniq_age
  }
  
  cohort <- as.integer(factor(year - age))
  Ip <- diag(1,ncx,ncx)
  
  Zc <- matrix(0, n0, ncoh) # cohort matrix
  Zc[cbind(1:n0, cohort)] <- 1

  
  ## 
  Xm <- matrix(0, n0, nobs*ncx)
  for(i in 1:nobs)
  {
    indexi <- indexy == uniq_index[i]
    Xm[indexi, (ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniq_index[i],]
  }
  
  ng <- length(unique(group.individual))
  W <- matrix(0, nobs, ng)
  W[cbind(1:nobs,group.individual)] <- 1
  W <- W %x% diag(1,ncx)
  
  ngc <- length(unique(group.cohort))
  Wc <- matrix(0, ncoh, ngc)
  Wc[cbind(1:ncoh, group.cohort)] <- 1
  Zc <- Zc%*%Wc
  
  Ux <- cbind(Xm%*%W, Zc) ## new X with group information
  
  ncu <- ncol(Ux)
  
  if(ng < nobs)
  {
    H1 <- matrix(0, nrow = 1, ncol = ncu) # constraint including beta 
    H1[,(ng*ncx+1):ncu] <- 1
    
    U1 <- matrix(0, ncu + 1,  ncol = ncu + 1) # linear model with constraint
    U1[1:ncu, 1:ncu] <- t(Ux)%*%Ux
    
    U1[ncu + 1, 1:ncu] <- H1
    U1[1:ncu, ncu+1] <- t(H1)
    thetaest <- (solve(U1)%*%c(t(Ux)%*%y,0))[1:ncu]
    
  }
  
  if(ng == nobs)
  {
    nz <- ncol(Zc)
    # three constraints
    H1 <- matrix(0, nrow = 3, ncol = ncu) # constraint including beta 
    H1[1,(ng*ncx+1):ncu] <- 1
    H1[2,(ng*ncx+1):ncu] <- 1:nz
    H1[3,(ng*ncx+1):ncu] <- (1:nz)^2
    
    U1 <- matrix(0, ncu + 3,  ncol = ncu + 3) # linear model with constraint
    U1[1:ncu, 1:ncu] <- t(Ux)%*%Ux
    
    U1[ncu + 1:3, 1:ncu] <- H1
    U1[1:ncu, ncu+1:3] <- t(H1)
    thetaest <- (solve(U1)%*%c(t(Ux)%*%y,0,0,0))[1:ncu]
  }
  
  betaest <- matrix(thetaest[1:(ng*ncx)],ncol = ncx, byrow=TRUE)
  etaest <- thetaest[-(1:(ng*ncx))]
  
  estimates <- Ux %*% thetaest
  
  sig2 <- sum((y - estimates)^2)/n0
  
  curve_est <- Xm%*%W%*%thetaest[1:(ng*ncx)]
  
  out <- list(betaest = betaest, etaest = etaest, estimates = estimates, curve = curve_est, sig2 = sig2)
  return(out)
}