#### a small simulation to see poential computational bias ##### 
### four groups, simulated from the existing model and data.... 
## ignore cohort effect ### sig2 = 0.08336409###
dat <- read.csv("/Users/wangx172/Dropbox/Tanja&XinW/Newdata/AggObese1990_2017.csv")
load("betaest4567.RData")
source("Gr_cohort/refit_cohort.R")

dat2 <- arrange(dat, AGE, IYEAR)
dat2 <- filter(dat2, AGE !=80)
year2 <- dat2$IYEAR
age2 <- dat2$AGE
nage2 <- length(unique(age2))
y2 <- dat2$PropObese 
x2 <- cbind(1, scale(dat2$IYEAR),scale((dat2$IYEAR - mean(dat2$IYEAR))^2))
ncoh <- nage2 + length(unique(year2)) - 1
ages = unique(dat2$AGE)

betaest40 = round(res_fitg4$betaest,3)
betaest40g = betaest40[group_coef4,]

sig2est = 0.083


bias_array4 = bias_array4r = array(0, dim = c(4,3,100))
bias_arrayall4 = bias_arrayall4r = array(0, dim = c(nage2, 3, 100))

relbias_array4 = relbias_array4r = array(0, dim = c(4,3,100))
relbias_arrayall4 = relbias_arrayall4r = array(0, dim = c(nage2, 3, 100))

for(mm in 1:100)
{
  set.seed(506 + mm)
  ysim = rep(0, length(y2))
  
  for(i in 1:nage2)
  {
    ysim[age2 == ages[i]] = x2[age2 == ages[i],] %*% betaest40[group_coef4[i],] 
  }
  
  ysim = ysim + rnorm(length(y2))*sqrt(sig2est)
  
  lam2c <- seq(0.05,5,by = 0.05)
  wmat <- matrix(0, nage2, nage2)
  for(i in 1:(nage2-1))
    for(j in (i+1):(nage2))
      wmat[i,j] <- abs(i - j)
  wmat <- wmat + t(wmat)
  ordervec <- wmat[lower.tri(wmat)]
  
  weights4 <- exp(0.75*(1-ordervec))
  betam04 <- cal_initialrx(indexy = age2,y = ysim,x = x2)
  
  bic_s4 = rep(0, length(lam2c)) 
  
  BICnt = function(obj, indexy, y, x)
  {
    betaest = obj$betaest
    groupest = obj$group
    ngroup = length(unique(groupest))
    n0 = length(y)
    
    uind = unique(indexy)
    meanest = rep(0,n0)
    for(i in 1:length(uind))
    {
      indi = indexy == uind[i]
      meanest[indi] = x[indi,] %*% betaest[groupest[i],]
    }
    
    BICvalue = log(sum((y - meanest)^2/n0)) + log(n0)/n0*(ngroup*ncol(x))
    return(BICvalue)
  }
  
  
  for(j in 1:length(lam2c))
  {
    res_s4j = Spgrrx(indexy = age2, y = ysim,x = x2,weights = weights4,
                     betam0 = betam04,lam = lam2c[j], maxiter = 1000)
    bic_s4[j] <-  BICnt(res_s4j,indexy = age2,y = ysim, x = x2)
  }
  
  res_fit4 = Spgrrx(indexy = age2, y = ysim,x = x2,weights = weights4,
                    betam0 = betam04,lam = lam2c[which.min(bic_s4)], maxiter = 1000)
  
  groupest = getgroup(res_fit4$deltam,n = nage2,tol = 0.005)
  betaest = res_fit4$betaest
  
  #### refit ##
  res_rfit4 = refit_group(year = dat2$IYEAR, age = age2, y = ysim, x = x2, group = groupest, model = "age")
  betaestr = matrix(res_rfit4$betaest, ncol = 3, byrow = TRUE)
  
  
  #### bias ###

  
  bias_array4[,,mm] = abs(betaest40 - betaestr)
  bias_array4r[,,mm] = abs(betaest40 - betaest)
  
  
  bias_arrayall4[,,mm] = abs(betaest40g - betaestr[groupest,])
  bias_arrayall4r[,,mm] = abs(betaest40g - betaest[groupest,])
  
  ### relative bias ###
  
  
  relbias_array4[,,mm]= abs(betaest40 - betaestr)/abs(betaest40)
  relbias_array4r[,,mm] = abs(betaest40 - betaest)/abs(betaest40)
  
  
  relbias_arrayall4[,,mm] =  abs(betaest40g - betaestr[groupest,])/abs(betaest40g)
  relbias_arrayall4r[,,mm] = abs(betaest40g - betaest[groupest,])/abs(betaest40g)
  
  print(mm)
  
}

par(mfrow = c(4,3))
for(i in 1:4)
{
  for(j in 1:3)
  {
    boxplot(bias_array4[i,j,],bias_array4r[i,j,])
  }
}

par(mfrow = c(4,3))
for(i in 1:4)
{
  for(j in 1:3)
  {
    boxplot(relbias_array4[i,j,],relbias_array4r[i,j,])
  }
}




