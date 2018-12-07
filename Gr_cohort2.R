####### Find some subgroups for both regression coefficients#######
#######and cohort effect for obesity data by Xin Wang ###### 
##### clustering regression coefficients and cohort effect together ####

library(Spgr)

######################### scad penalty ###################
sfun <- function(x, th)
{
  xn <- sqrt(sum(x^2))
  thval <- 1 - th/xn
  thval*((thval) >0)*x
}

# the solution of the scad penalty 
scad <- function(x,lam,nu,gam)
{
  temp1 <- lam/nu
  temp2 <- gam * lam
  xn <- sqrt(sum(x^2))
  
  if(xn <= lam + temp1)
  {
    z <- sfun(x, temp1)
  }else if(xn <= temp2 & xn >= lam + temp1)
  {
    z <- sfun(x, temp2/((gam-1)*nu))/(1- 1/((gam - 1 )*nu))
  }else{
    z <- x
  }
  
  return(z)
  
}
###########################################################


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

# the data year and age should sort increasingly.


Gr_cohort2 <- function(year, age, y, x, betam0, model = "year", weights,
                      lam1 = 0.1, lam2 = 0.1, nu = 1, gam = 3, 
                      maxiter = 1000, tolabs = 1e-4, tolrel = 1e-2)
{
  n0 <- length(y) # number of total observation
  ncx <- ncol(x) # dimension of x
  
  uniq_year <- unique(year) # unique index sort increasingly 
  uniq_age <- unique(age) # unique age sort increasingly 
  
  nyear <- length(uniq_year)
  nage <- length(uniq_age)
  ncoh <- nyear + nage - 1 # number of cohort effects
  
  # cluster year curves
  if(model == "year")
  {
    nobs <- nyear ## number of individuals
    nrep <- nage # number of observations for each indiviual
    indexy <- year
    uniq_index <- uniq_year
    
  }
  
  # cluster age curves 
  if(model == "age")
  {
    nobs <- nage
    nrep <- nyear
    indexy <- age
    uniq_index <- uniq_age
  }
  
  Ip <- diag(1,ncx,ncx)
  
  npair <- nobs * (nobs - 1)/2 # number of pairs for regression coefficients
  npairc <- ncoh * (ncoh - 1)/2 # number of pairs for cohort effect 
  
  
  # Differene matrix for coefficients
  D <- matrix(0,npair,nobs)  # ei - ej
  for(j in 1:(nobs-1))
  {
    indexj <- (nobs-1 + nobs-j+1)*(j-1)/2
    indexvj <- indexj + (1:(nobs-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nobs)] <- -1
  }
  
  ## Difference matrix for cohort effect
  Dc <- matrix(0, npairc, ncoh)
  for(j in 1:(ncoh-1))
  {
    indexj <- (ncoh-1 + ncoh-j+1)*(j-1)/2
    indexvj <- indexj + (1:(ncoh-j))
    Dc[indexvj,j] <- 1
    Dc[cbind(indexvj,(j+1):ncoh)] <- -1
  }
  
  
  AtA <- t(D)%*%D %x% Ip
  
  cohort <- as.integer(factor(year - age)) # cohort effect 
  Zc <- matrix(0, n0, ncoh) # cohort matrix
  Zc[cbind(1:n0, cohort)] <- 1
  
  Hm <- matrix(0, 1, ncoh) # constraints matrix
  Hm[1,] <- rep(1,ncoh)

  
  #### transformation including W matrix 
  Xm <- X <- matrix(0, n0, nobs*ncx)
  Zcm <- Zc
  ym <- rep(0, n0)
  for(i in 1:nobs)
  {
    indexi <- indexy == uniq_index[i]
    Xm[indexi, (ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniq_index[i],]/sqrt(nrep)
    X[indexi, (ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniq_index[i],]
    ym[indexi] <- y[indexi]/sqrt(nrep)
    Zcm[indexi,] <- Zc[indexi,]/sqrt(nrep)
  }
  
  
  XtX_inv <- solve(t(Xm) %*% Xm + nu * AtA)
  ZtZ_inv <- solve(t(Zcm) %*% Zcm + nu * t(Hm) %*% Hm + nu * t(Dc)%*% Dc)
  Xty <- t(Xm) %*%ym
  Zty <- t(Zcm) %*%ym
  XtZ <- t(Xm) %*% Zcm
  ZtX <- t(XtZ)
  
  reg_b1 <- XtX_inv %*% Xty
  reg_eta1 <- ZtZ_inv %*% Zty
  
  
  ##### some initials ####
  betam <- betam0
  betanew <- c(t(betam)) # vector form 
  
  ## initial for eta ###
  
  etanew <-  solve(t(Zcm)%*%Zcm) %*% t(Zcm)%*%(ym - Xm%*%betanew)
  etanew <- etanew - mean(etanew)
  
  
  deltam.old <- t(D %*% betam)
  alpm.old <-  Dc%*% etanew
  
  vm  <-  matrix(0, ncx, nobs*(nobs-1)/2)
  vmc <- rep(0, npairc)
  vh <-  0 
  
  flag <- 0
  
  for(m in 1:maxiter)
  {
    # update beta and eta 
    betanew <- reg_b1 - XtX_inv %*% XtZ %*% etanew + nu* XtX_inv %*% c((deltam.old -  vm/nu) %*% D)
    etanew <- reg_eta1 + ZtZ_inv %*% ( - ZtX %*% betanew - t(Hm) %*% vh + 
                                         nu * t(Dc) %*% (alpm.old - vmc/nu))
    
    betam <- matrix(betanew, nobs, ncx, byrow = TRUE) # matrix 
    betadiff <- t(D %*% betam)
    psim <- betadiff + vm/nu
    
    etadiff <- Dc %*% etanew
    taum <- etadiff + vmc/nu
    
    
    # update delta and alpm
    deltam <- sapply(1:ncol(psim),function(xx) scad(psim[,xx],weights[xx]*lam1 ,nu,gam))
    alpm <- sapply(1:length(taum),function(xx) scad(taum[xx], lam2 ,nu,gam))
    
    # update v
    Abd <- betadiff - deltam
    Heta <- Hm %*% etanew
    Dea <- etadiff - alpm
    
    vm <- vm + nu * Abd
    vmc <- vmc + nu * Dea
    vh <- vh + nu * Heta
    
    # convergence
    rm <- sqrt(sum(Abd^2) + Heta^2)
    sm <- nu*sqrt(sum(((deltam - deltam.old)%*%D)^2) + sum((t(Dc)%*%(alpm - alpm.old))^2))
    
    tolpri <- tolabs*sqrt(npair*ncx + npairc + 1) + 
      tolrel*max(sqrt(sum(betadiff^2) + Heta^2 + sum(etadiff^2)),sqrt(sum(deltam^2) + sum(alpm^2)))
    toldual <- tolabs*sqrt(nobs*ncx + ncoh) + tolrel* sqrt(sum((vm %*% D)^2) + vh^2*ncoh + sum((t(Dc)%*% vmc)^2))
    
    deltam.old <- deltam
    alpm.old <- alpm
    
    if(rm <= tolpri & sm <= toldual){break}
  }
  
  
  if(m == maxiter) {flag <- 1}
  
  groupest <- getgroup(deltam = deltam,n = nobs) 
  # getgroup is a function Spgr to find the group information based on estimated delta
  ngest <- length(unique(groupest))
  
  ## group information for cohort 
  groupestc <- getgroup(deltam = matrix(alpm,nrow = 1), n = ncoh)
  ngestc <- length(unique(groupestc))
  # group coefficients
  if(ncx ==1){
    alpest <-  matrix(by(betam, groupest, colMeans,simplify = TRUE), nrow= 1)
  }else{  alpest <- do.call("rbind",by(betam, groupest, colMeans,simplify = TRUE))}
  
  betaest <- alpest[groupest,]
  
  betac <- c(t(betaest))
  BICvalue <-  log(sum(ym - Xm %*% betac - Zcm %*% etanew)^2/nobs) + log(nobs)/nobs*(ngest * ncx + ngestc)
  
  BICvalue2 <- log(sum(y - X%*% betac - Zc%*% etanew)^2/n0) + log(n0)/n0*(ngest * ncx + ngestc)
  
  BICc <- log(sum(ym - Xm %*% betac - Zcm %*% etanew)^2/nobs) + log(log(nobs*ncx))*log(nobs)/nobs*(ngest * ncx + ngestc)
  
  BICc2 <- log(sum(y - X%*% betac - Zc%*% etanew)^2/n0) + log(log(nobs*ncx + ncoh))*log(n0)/n0*(ngest * ncx + ngestc)
  
  BIClog <- log(sum(y - X%*% betac - Zc%*% etanew)^2/n0) + log(nobs*ncx + ncoh)*log(n0)/n0*(ngest * ncx + ngestc)
  
  resi <- sum(y - X%*% betac - Zc%*% etanew)^2/n0
  
  outls <- list(betaest = betaest, etaest = etanew, betam = betam, alpest = alpest,
                group = groupest, groupc = groupestc,
                deltam = deltam,  alpm = alpm, 
                BIC = BICvalue, BIC2 = BICvalue2,BICc = BICc, BICc2 = BICc2, BIClog = BIClog, resi = resi,
                rm = rm, sm = sm, tolpri = tolpri, toldual = toldual,
                flag = flag, niteration = m)
  
  return(outls)
}





