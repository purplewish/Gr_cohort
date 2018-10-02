####### Find some subgroups for obesity data by Xin Wang####

# y is the proportion
# age is the age 
# x is the covariate, which are based on ages
# year is "year" 
# weights: weight for each pair
# betam0: initial value
# lam: tuning parameter
# maxiter: maximum number of iterations
# tolabs and tolrel are two tolerance criteria in ADMM
# model: model can be year or age. If year is specified, clustering is for year

dat <- read.csv("Rfiles/CBD-O/FullObeseYear.csv")
year <- dat$IYEAR
age <- dat$AGE
y <- dat$PropObese
x <- cbind(1, scale(dat$AGE), scale(dat$AGE^2))
nu <- 1
gam <- 3

group <- rep(1:2,c(12,13))

temp <- refit_cohort(year = year,age = age,x = x,group = group)
Xm <- temp$Xm


Gr_cohort <- function(year, age, y, x, betam0, model = "year", weights,
                      lam = 0.5, nu = 1, gam = 3, lam = 0.5,
                      maxiter = 1000, tolabs = 1e-4, tolrel = 1e-2)
{
  n0 <- length(y) # number of total observation
  ncx <- ncol(x)
  
  uniq_year <- unique(year) # unique index sort increasingly 
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
  
  D <- matrix(0,nobs*(nobs-1)/2,nobs)  # ei - ej
  for(j in 1:(nobs-1))
  {
    indexj <- (nobs-1 + nobs-j+1)*(j-1)/2
    indexvj <- indexj + (1:(nobs-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nobs)] <- -1
  }
  
  AtA <- t(D)%*%D %x% Ip
  
  
  Zc <- matrix(0, n0, ncoh) # cohort matrix
  Zc[cbind(1:n0, cohort)] <- 1
  
  Hm <- matrix(0, 3, ncoh) # constraints matrix
  Hm[1,] <- rep(1,ncoh)
  Hm[2,]<-  1:ncoh
  Hm[3,] <- (1:ncoh)^2
  
  
  
  #### transformation, reoder the data based on uniq_index
  Xm <- matrix(0, n0, nobs*ncx)
  Zcm <- Zc
  ym <- rep(0, n0)
  for(i in 1:nobs)
  {
    indexi <- indexy == uniq_index[i]
    Xm[indexi, (ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniq_index[i],]/sqrt(nrep)
    ym[indexi] <- y[indexi]/sqrt(nrep)
    Zcm[indexi,] <- Zc[indexi,]/sqrt(nrep)
  }
  
  
  XtX_inv <- solve(t(Xm) %*% Xm + nu * AtA)
  ZtZ_inv <- solve(t(Zcm) %*% Zcm + nu * t(Hm) %*% Hm)
  Xty <- t(Xm) %*%ym
  Zty <- t(Zcm) %*%ym
  XtZ <- t(Xm) %*% Zcm
  ZtX <- t(XtZ)
  
  reg_b1 <- XtX_inv %*% Xty
  reg_eta1 <- ZtZ_inv %*% Zty
  ##### some initials ####
  deltam.old <- t(D %*% betam0)
  betam <- betam0
  betanew <- c(t(betam))
  eta <- rep(0,ncoh)
  
  vm  <-  matrix(0, ncx, nobs*(nobs-1)/2)
  vh <- rep(0, 3)
  
  flag <- 0
  
  for(m in 1:maxiter)
  {
    betanew <- reg_b1 - XtX_inv %*% XtZ %*% eta + nu* XtX_inv %*% c((deltam.old -  vm/nu) %*% D)
    etanew <- reg_eta1 - ZtZ_inv %*% ZtX %*% betanew - t(Hm) %*% vh
    
    betadiff <- t(D %*% betam)
    psim <- betadiff + vm/nu
    deltam <- sapply(1:ncol(psim),function(xx) scad(psim[,xx],cvec[xx]*lam,nu,gam))
    vm <- vm + nu * (betadiff - deltam)
    vh <- vh + nu * Hm %*% etanew 
  }
  
  
  
  
}


