#### use simulated data #####

### two groups ##### 
source("Gr_cohort/refit_cohort.R")

beta20 = matrix(c(-0.5,-0.5,0.5,0.5), ncol = 2, byrow = TRUE)
sig2 = 0.4
n = 60 ### number of subjects 
nrep = 30 ### number of replicates for each subject
group2 = rep(1:2, each = 30)
beta20g = beta20[group2,]

nsim = 200

bias_array2 = bias_array2r = array(0, dim = c(2,2,nsim))
bias_arrayall2 = bias_arrayall2r = array(0, dim = c(n, 2, nsim))

relbias_array2 = relbias_array2r = array(0, dim = c(2,2,nsim))
relbias_arrayall2 = relbias_arrayall2r = array(0, dim = c(n, 2, nsim))

indexy2 = rep(1:n, each = nrep)

for(mm in 1:nsim)
{
  set.seed(506 + mm)
  ysim = rep(0, n*nrep)
  
  x2 = matrix(rnorm(n*nrep*2), ncol = 2)
  
  for(i in 1:n)
  {
    ysim[indexy2 == i] = x2[indexy2 == i,] %*% beta20[group2[i],] 
  }
  
  ysim = ysim + rnorm(n*nrep)*sqrt(sig2)
  
  lam2c <- seq(0.05,1,by = 0.01)
  
  betam04 <- cal_initialrx(indexy = indexy2,y = ysim,x = x2)
  
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
  
  wts = rep(1, n*(n-1)/2)
  for(j in 1:length(lam2c))
  {
    res_s4j = Spgrrx(indexy = indexy2, y = ysim,x = x2,weights = wts,
                     betam0 = betam04,lam = lam2c[j], maxiter = 1000)
    bic_s4[j] <-  BICnt(res_s4j,indexy = indexy2,y = ysim, x = x2)
  }
  
  res_fit4 = Spgrrx(indexy = indexy2, y = ysim,x = x2,weights = wts,
                    betam0 = betam04,lam = lam2c[which.min(bic_s4)], maxiter = 1000)
  
  groupest = getgroup(res_fit4$deltam,n = n,tol = 0.005)
  betaest = res_fit4$betaest
  
  #### refit ##
  res_rfit4 = refit_group(year = indexy2, age = indexy2, y = ysim, x = x2, group = groupest, model = "age")
  betaestr = matrix(res_rfit4$betaest, ncol = 2, byrow = TRUE)
  
  
  #### bias ###
  
  if(length(unique(groupest))==2)
  {
    bias_array2[,,mm] = abs(beta20 - betaestr)
    bias_array2r[,,mm] = abs(beta20 - betaest)
    
    
    bias_arrayall2[,,mm] = abs(beta20g - betaestr[groupest,])
    bias_arrayall2r[,,mm] = abs(beta20g - betaest[groupest,])
    
    ### relative bias ###
    
    
    relbias_array2[,,mm]= abs(beta20 - betaestr)/abs(beta20)
    relbias_array2r[,,mm] = abs(beta20 - betaest)/abs(beta20)
    
    
    relbias_arrayall2[,,mm] =  abs(beta20g - betaestr[groupest,])/abs(beta20g)
    relbias_arrayall2r[,,mm] = abs(beta20g - betaest[groupest,])/abs(beta20g)
    
  }
  
 
  print(mm)
  
}

save(bias_array2,bias_array2r,bias_arrayall2,bias_arrayall2r,
     relbias_array2,relbias_array2r,relbias_arrayall2, relbias_arrayall2r, file = "sim_bias.RData")


load("sim_bias.RData")


boxplot(apply(bias_array2,3,sum),apply(bias_array2r,3,sum))

par(mfrow = c(2,2))
for(i in 1:2)
{
  for(j in 1:2)
  {
    boxplot(bias_array2[i,j,],bias_array2r[i,j,])
  }
}

par(mfrow = c(2,2))
for(i in 1:2)
{
  for(j in 1:2)
  {
    boxplot(relbias_array2[i,j,],relbias_array2r[i,j,])
  }
}




