####### Obesity data ####

# y is the proportion
# age is the age 
# x is the covariate, which are based on ages
# year is "year" 
# betam0: initial value
# lam: tuning parameter
# maxiter: maximum number of iterations
# tolabs and tolrel are two tolerance criteria in ADMM

Gr_cohort <- function(year, age, y, x, betam0, model = "year",
                      lam = 0.5, nu = 1, gam = 3, lam = 0.5,
                      maxiter = 1000, tolabs = 1e-4, tolrel = 1e-2)
{
  n0 <- length(y) # number of total observation
  ncx <- ncol(x)
  
  uniq_year <- unique(year)
  uniq_age <- unique(age)
  
  nyear <- length(uniq_year)
  nage <- length(uniq_age)
  
  nobs <- length(uniqxy)
  npair <- nobs*(nobs-1)/2
  
  Ip <- diag(1,ncx,ncx)
  
  D <- matrix(0,nobs*(nobs-1)/2,nobs)
  for(j in 1:(nobs-1))
  {
    indexj <- (nobs-1 + nobs-j+1)*(j-1)/2
    indexvj <- indexj + (1:(nobs-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nobs)] <- -1
  }
  
  AtA <- t(D)%*%D %x% Ip
}