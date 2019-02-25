###### Gr_cohort_only: Find subgroups for cohort given the group information of regression coefficients ####
####### Gr_coef_only: Find subgroups for regression coefficients given the group information of cohort #######

source("Gr_cohort/scad.R")

Gr_cohort_only <- function(year, age, y, x, model = "year",group.individual,
                                          lam2 = 0.1, nu = 1, gam = 3, 
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
  
  npairc <- ncoh * (ncoh - 1)/2 # number of pairs for cohort effect 
  
  ## Difference matrix for cohort effect
  Dc <- matrix(0, npairc, ncoh)
  for(j in 1:(ncoh-1))
  {
    indexj <- (ncoh-1 + ncoh-j+1)*(j-1)/2
    indexvj <- indexj + (1:(ncoh-j))
    Dc[indexvj,j] <- 1
    Dc[cbind(indexvj,(j+1):ncoh)] <- -1
  }
  
  
  cohort <- as.integer(factor(year - age)) # cohort effect 
  Zc <- matrix(0, n0, ncoh) # cohort matrix
  Zc[cbind(1:n0, cohort)] <- 1
  
  
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
  
  ng <- length(unique(group.individual))
  
  if(ng < nobs)
  {
    W <- matrix(0, nobs, ng)
    W[cbind(1:nobs,group.individual)] <- 1
    W <- W %x% diag(1,ncx)
    Xm <- Xm%*%W 
    X <- X%*%W
    
    Hm <- matrix(0, 1, ncoh) # constraints matrix
    Hm[1,] <- rep(1,ncoh)
  }else{
    Hm <- matrix(0, 3, ncoh) # constraints matrix
    Hm[1,] <- rep(1,ncoh)
    Hm[2,]<-  ncoh*(1:ncoh)/sum(1:ncoh)
    Hm[3,] <- ncoh*(1:ncoh)^2/(sum((1:ncoh)^2))
  }

  nh <- nrow(Hm)

  XtX_inv <- solve(t(Xm)%*%Xm)
  
  Qx <- diag(n0) - Xm %*% XtX_inv%*%t(Xm)
  
  ZtZ_inv <- solve(t(Zcm) %*% Qx%*%Zcm + nu * t(Hm) %*% Hm + nu * t(Dc)%*% Dc)
  Zty <- t(Zcm) %*%Qx%*%ym
  reg_eta1 <- ZtZ_inv %*% Zty
  
  ##### some initials ####
  betanew <- XtX_inv %*%t(Xm)%*%ym   ##vector form

  ## initial for eta ###
  
  etanew <-  solve(t(Zcm)%*%Zcm) %*% t(Zcm)%*%(ym - Xm%*%betanew)
  etanew <- etanew - mean(etanew)
  
  alpm.old <-  Dc%*% etanew
  
  vmc <- rep(0, npairc)
  vh <-  rep(0, nh)
  
  flag <- 0
  
  
  for(m in 1:maxiter)
  {
    # update eta 
    etanew <- reg_eta1 + ZtZ_inv %*% (- t(Hm) %*% vh + 
                                         nu * t(Dc) %*% (alpm.old - vmc/nu))
    
    etadiff <- Dc %*% etanew
    taum <- etadiff + vmc/nu
    
    # alpm
    alpm <- sapply(1:length(taum),function(xx) scad(taum[xx], lam2 ,nu,gam))
    
    # update v
    Heta <- Hm %*% etanew
    Dea <- etadiff - alpm
    
    vmc <- vmc + nu * Dea
    vh <- vh + nu * Heta
    
    # convergence
    rm <- sqrt(sum(Dea)^2 + sum(Heta)^2)
    sm <- nu*sqrt(sum((t(Dc)%*%(alpm - alpm.old))^2))
    
    tolpri <- tolabs*sqrt(npairc + nh) + 
      tolrel*max(sqrt(sum(Heta)^2 + sum(etadiff^2)),sqrt(sum(alpm^2)))
    toldual <- tolabs*sqrt(ncoh) + tolrel* sqrt(sum((t(Hm)%*%vh + t(Dc)%*% vmc)^2))
  
    alpm.old <- alpm
    
    if(rm <= tolpri & sm <= toldual){break}
  }
  
  if(m == maxiter) {flag <- 1}

  betanew <- XtX_inv %*%t(Xm)%*%(ym - Zcm%*%etanew)  
  
  betam <- matrix(betanew, ng, ncx, byrow = TRUE)
  
  
  ## group information for cohort 
  groupestc <- getgroup(deltam = matrix(alpm,nrow = 1), n = ncoh,tol = 1e-8)
  ngestc <- length(unique(groupestc))

  
  BICvalue2 <- log(sum(y - X%*% betanew - Zc%*% etanew)^2/n0) + log(n0)/n0*(ng*ncx + ngestc )
  
  BICc2 <- log(sum(y - X%*% betanew - Zc%*% etanew)^2/n0)  + log(log(ng*ncx + ncoh))*log(n0)/n0*(ng*ncx + ngestc)
  
  outls <- list(etaest = etanew, betam = betam,  groupc = groupestc, alpm = alpm, BIC2 = BICvalue2, BICc2 = BICc2,
                rm = rm, sm = sm, tolpri = tolpri, toldual = toldual,
                flag = flag, niteration = m)
  return(outls)
}


##### given group information of cohort update beta 
#Gr_coef_only 

Gr_coef_only <- function(year, age, y, x, model = "year", group.cohort, betam0,
                         ws, lam = 0.1, nu = 1, gam = 3, 
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
  
  npair <- nobs*(nobs - 1)/2
  
  D <- matrix(0,npair,nobs)  # ei - ej
  for(j in 1:(nobs-1))
  {
    indexj <- (nobs-1 + nobs-j+1)*(j-1)/2
    indexvj <- indexj + (1:(nobs-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nobs)] <- -1
  }
  
  AtA <- t(D)%*%D %x% Ip
  
  cohort <- as.integer(factor(year - age)) # cohort effect 
  Zc <- matrix(0, n0, ncoh) # cohort matrix
  Zc[cbind(1:n0, cohort)] <- 1
  
  
  
  Xm <- X <- matrix(0, n0, nobs*ncx)
  Zcm <- Zc
  ym <- rep(0, n0)
  for(i in 1:nobs)
  {
    indexi <- indexy == uniq_index[i]
    Xm[indexi, (ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniq_index[i],]/sqrt(nrep)
    ym[indexi] <- y[indexi]/sqrt(nrep)
    Zcm[indexi,] <- Zc[indexi,]/sqrt(nrep)
    X[indexi, (ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniq_index[i],]
  }
  
  ngc <- length(unique(group.cohort))
  Wc <- matrix(0, ncoh, ngc)
  Wc[cbind(1:ncoh, group.cohort)] <- 1
  Zc <- Zc%*%Wc
  Zcm <- Zcm %*% Wc
  
  Hm <- matrix(0, 1, ngc) # constraints matrix
  Hm[1,] <- rep(1,ngc)
  
  
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
  betanew <- c(t(betam)) # vector form 
  etanew <- rep(0,ngc)
  
  vm  <-  matrix(0, ncx, npair)
  vh <-  0 
  
  flag <- 0
  
  for(m in 1:maxiter)
  {
    # update beta and eta 
    betanew <- reg_b1 - XtX_inv %*% XtZ %*% etanew + nu* XtX_inv %*% c((deltam.old -  vm/nu) %*% D)
    etanew <- reg_eta1 - ZtZ_inv %*% ZtX %*% betanew - ZtZ_inv %*%t(Hm) %*% vh
    
    betam <- matrix(betanew, nobs, ncx, byrow = TRUE) # matrix 
    betadiff <- t(D %*% betam)
    psim <- betadiff + vm/nu
    
    # update delta 
    deltam <- sapply(1:ncol(psim),function(xx) scad(psim[,xx],ws[xx]*lam,nu,gam))
    
    # update v
    Abd <- betadiff - deltam
    Heta <- Hm %*% etanew
    
    vm <- vm + nu * Abd
    vh <- vh + nu * Heta
    
    # convergence
    rm <- sqrt(sum(Abd^2) + Heta^2)
    sm <- nu*sqrt(sum(((deltam - deltam.old)%*%D)^2))
    
    tolpri <- tolabs*sqrt(npair*ncx + 1) + tolrel*max(sqrt(sum(betadiff^2) + Heta^2),sqrt(sum(deltam^2)))
    toldual <- tolabs*sqrt(nobs*ncx + ngc) + tolrel* sqrt(sum((vm %*% D)^2) + vh^2*ngc)
    
    deltam.old <- deltam
    
    if(rm <= tolpri & sm <= toldual){break}
  }
  
  if(m == maxiter) {flag <- 1}
  
  groupest <- getgroup(deltam = deltam,n = nobs) 
  # getgroup is a function Spgr to find the group information based on estimated delta
  ngest <- length(unique(groupest))
  
  # group coefficients
  if(ncx ==1){
    alpest <-  matrix(by(betam, groupest, colMeans,simplify = TRUE), nrow= 1)
  }else{  alpest <- do.call("rbind",by(betam, groupest, colMeans,simplify = TRUE))}
  
  betaest <- alpest[groupest,]
  
  betac <- c(t(betaest))

  BICvalue2 <- log(sum(y - X%*%betac - Zc%*% etanew)^2/n0) + log(n0)/n0*(ngest*ncx + ngc)
  
  BICc2 <- log(sum(y - X%*%betac- Zc%*% etanew)^2/n0)  + 0.2*log(log(nobs*ncx + ngc))*log(n0)/n0*(ngest*ncx + ngc)
  
  outls <- list(betaest = betaest, etaest = etanew, betam = betam, alpest = alpest,
                group = groupest, deltam = deltam, BIC2 = BICvalue2, BICc2 = BICc2,
                rm = rm, sm = sm, tolpri = tolpri, toldual = toldual,
                flag = flag, niteration = m)
  return(outls)
}