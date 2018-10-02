############# refit algorithm by Xin Wang#######
# refit based on a given group information
# index should be ordered increasingly
# group is the predefined group information, the order is the same as index (year or age)
# model can be year or age. If year is specified, group is for year

refit_cohort <- function(year, age, x, group, model = "year")
{
  
  n0 <- length(y) # number of total observation
  ncx <- ncol(x)
  
  uniq_year <- unique(year) 
  uniq_age <- unique(age)
  
  nyear <- length(uniq_year)
  nage <- length(uniq_age)
  ncoh <- nyear + nage - 1 # number of cohort effects
  
  cohort <- as.integer(factor(year - age))
  
  Ip <- diag(1,ncx,ncx)
  
  if(model == "year")
  {
    nobs <- nyear ## number of individuals
    nrep <- nage # number of observations for each indiviual
    indexy <- year
    uniq_index <- uniq_year
  }
 
  Xm <- matrix(0, n0, nobs*ncx)
  ym <- rep(0, n0)
  
  for(i in 1:nobs)
  {
    Xm[indexy == uniq_index[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniq_index[i],]/sqrt(nrep)
    ym[indexy == uniq_index[i]] <- y[indexy == uniq_index[i]]/sqrt(nrep)
  }
  
  ng <- length(unique(group))
  W <- matrix(0, nobs, ng)
  W[cbind(1:nobs,group)] <- 1
  W <- W %x% diag(1,ncx)
  Ux <- Xm%*%W
  est <- solve(t(Ux)%*%Ux)%*%t(Ux)%*%ym
  
  sig2 <- sum((ym - Ux%*%est)^2)/nobs
  
  betaest <- matrix(est, ng, ncx, byrow = TRUE)
  
  #out <- list(est = betaest, sig2 = sig2)
  out <- list(est = betaest, sig2 = sig2, Xm = Xm)
  return(out)
}


#refit_cohort(year = year,age = age,x = x,group = group)
