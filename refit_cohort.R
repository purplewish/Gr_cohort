############# refit algorithm by Xin Wang#######
# refit based on a given group information
# year and age should be ordered increasingly
# group is the predefined group information, the order is the same as index (year or age)
# model can be year or age. If year is specified, group is for year
# since equal size the weighted square is not necessary 



refit_cohort <- function(year, age, x, group, model = "year")
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
    index < age
    uniq_index <- uniq_age
  }
  
  cohort <- as.integer(factor(year - age))
  Ip <- diag(1,ncx,ncx)
  
  Zc <- matrix(0, n0, ncoh) # cohort matrix
  Zc[cbind(1:n0, cohort)] <- 1
  
  Hm <- matrix(0, 1, ncoh) # constraints matrix
  Hm[1,] <- rep(1,ncoh)
  
  
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
  H1 <- matrix(0, nrow = 1, ncol = ncol(Ux)) # constraint including beta 
  H1[,(ng*ncx+1):ncol(Ux)] <- rnorm(103-3*6)
  
  U_inv <- solve(t(Ux)%*%Ux + t(H1)%*%H1)
  
  fit1 <- lm(c(y,0)~ 0+  rbind(Ux,H1))
  
  thetaest <- (U_inv - U_inv%*%t(H1)%*% solve(H1%*%U_inv%*%t(H1))%*%H1%*%U_inv)%*%t(Ux)%*%y # solution 

  
  betaest <- matrix(thetaest[1:(ng*ncx)],ncol = ncx, byrow=TRUE)
  etaest <- thetaest[-(1:(ng*ncx))]

  sig2 <- sum((ym - Ux%*%thetaest)^2)/nobs

  out <- list(betaest = betaest, etaest = etaest, sig2 = sig2)
  return(out)
}


