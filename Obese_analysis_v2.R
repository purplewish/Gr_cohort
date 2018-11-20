####### fit coefficients and cohort effect together ####

library(Spgr)
library(plyr)
library(dplyr)
library(ggplot2)
setwd("Research/Obesity/")
source("Gr_cohort/Gr_cohort.R")
source("Gr_cohort/refit_cohort.R")
source("Gr_cohort/Gr_cohort2.R")
dat <- read.csv("/Users/wangx172/Dropbox/Tanja&XinW/Rfiles/CBD-O/FullObeseYear1.csv")

#source("/Users/wangx172/Research/Obesity/Gr_cohort/Gr_cohort.R")
#source("/Users/wangx172/Research/Obesity/Gr_cohort/refit_cohort.R")

#source("C:/Users/miljkot/Dropbox/Tanja&XinW/Rfiles/Gr_cohort/Gr_cohort.R")
#source("C:/Users/miljkot/Dropbox/Tanja&XinW/Rfiles/Gr_cohort/refit_cohort.R")

#source("/Users/wangx172/Research/Obesity/Gr_cohort/Gr_cohort.R")
#source("/Users/wangx172/Research/Obesity/Gr_cohort/refit_cohort.R")
#dat <- read.csv("C:/Users/miljkot/Dropbox/Tanja&XinW/Rfiles/CBD-O/FullObeseYear1.csv")


###### find group for years #####
dat <- arrange(dat, IYEAR, AGE)
year <- dat$IYEAR
nyear <- length(unique(year))
age <- dat$AGE
uage <- unique(dat$AGE)
uyear <- unique(dat$IYEAR)

ncoh <- nyear + length(unique(age)) - 1

x <- cbind(1, scale(dat$AGE), scale(dat$AGE^2))
nu <- 1
gam <- 3
weights <- rep(1, nyear*(nyear-1))
betam01 <- cal_initialrx(indexy = dat$IYEAR,y = dat$PropObese,x = x) # from Spgr package 



lam1 <- seq(0.001,0.03,by = 0.001)
lam2 <- seq(0.001,0.004,by = 0.0005)

bic1 <- rep(0, length(lam2))
beta_array1 <- array(0, dim = c(nyear,3,length(lam2)))
groupmat1 <- matrix(0, nyear, length(lam2))

etamat1 <- groupeta1 <- matrix(0, ncoh, length(lam2))


for(j2 in 1:length(lam2))
{
  betam1 <- betam01
  
  bic1j <- rep(0, length(lam1))
  beta_arrayj <- array(0, dim = c(nyear,ncol(x),length(lam1)))
  groupmatj <- matrix(0, nyear, length(lam1))
 
  etamatj <- matrix(0, ncoh, length(lam1))
  groupetaj  <- matrix(0, ncoh, length(lam1))
  
  for(j1 in 1:length(lam1))
  {
    resj <- Gr_cohort2(year = year,age = age,y = dat$PropObese,x = x,betam0 = betam1,model = "year",weights = weights,lam1 = lam1[j1],lam2 = lam2[j2])
    betam1 <- resj$betaest
    bic1j[j1] <- resj$BIC2
    beta_arrayj[,,j1] <- resj$betaest
    groupmatj[,j1] <- resj$group
    
    etamatj[,j1] <- resj$etaest
    groupetaj[,j1] <- resj$groupc
  }
  
  indj <- which.min(bic1j)
  bic1[j2] <- bic1j[indj]
  
  beta_array1[,,j2] <- beta_arrayj[,,indj]
  groupmat1[,j2] <- groupmatj[,indj]
  etamat1[,j2] <- etamatj[,indj]
  groupeta1[,j2] <- groupetaj[,j2]
  
  print(j2)
}

plot(beta_array1[,,3])
plot(etamat1[,3])

groupmat1[,3]
groupeta1[,3]

xm <- dat$AGE - mean(uage)
sigmax <- sum((uage - mean(uage))^2)/length(uage)

x0 <- cbind(1, xm, xm^2 - sigmax)

res_c1 <- refit_cohort2(year = year, age = age,y = dat$PropObese, x = x0,group.individual = groupmat1[,1],group.cohort = groupeta1[,3],model = "year")


####### for age ######
### weights are considered 

dat2 <- arrange(dat, AGE, IYEAR)
year2 <- dat2$IYEAR
age2 <- dat2$AGE
nage2 <- length(unique(age2))
y2 <- dat2$PropObese 
x2 <- cbind(1, scale(dat2$IYEAR),scale(dat2$IYEAR^2))
ncoh <- nage2 + length(unique(year2)) - 1

# weights for age is defined based on the age difference 
wmat <- matrix(0, nage2, nage2)
for(i in 1:(nage2-1))
  for(j in (i+1):(nage2))
    wmat[i,j] <- abs(i - j)
wmat <- wmat + t(wmat)
ordervec <- wmat[lower.tri(wmat)]

betam02 <- cal_initialrx(indexy = age2,y = y2,x = x2) # from Spgr package 


lam1 <- seq(0.01,0.4,by = 0.01)
lam2 <- seq(0.01,0.1,by = 0.02)
alp <- c(0.25,0.5,0.75,1,1.25,1.5)

bic2 <- rep(0, length(alp))
beta_array2 <- array(0, dim = c(nage2,3,length(alp)))
groupmat2 <- matrix(0, nage2, length(alp))
etamat2 <- groupeta2 <- matrix(0, ncoh, length(alp))


for(k in 1:length(alp))
{
  weightsk <- exp(alp[k]*(1-ordervec))
 
  bic2k <- rep(0, length(lam2))
  
  beta_arrayk2 <- array(0, dim = c(nage2,ncol(x2),length(lam2)))
  groupmatk2 <- matrix(0, nage2, length(lam2))
  
  etamatk2 <- matrix(0, ncoh, length(lam2))
  groupetak2  <- matrix(0, ncoh, length(lam2))
  
  for(j2 in 1:length(lam2))
  {
    bic2j <- rep(0, length(lam1))

    beta_arrayj2 <- array(0, dim = c(nage2,ncol(x2),length(lam1)))
    groupmatj2 <- matrix(0, nage2, length(lam1))
    
    etamatj2 <- matrix(0, ncoh, length(lam1))
    groupetaj2  <- matrix(0, ncoh, length(lam1))
    
    resi2 <- rep(0, length(lam1))
    for(j1 in 1:length(lam1))
    {
      resj <- Gr_cohort2(year = year2,age = age2,y = scale(dat$PropObese),x = x2,
                         betam0 = betam02,model = "age",weights = weightsk,
                         lam1 = lam1[j1],lam2 = lam2[j2],maxiter = 2000)
      
      resi2[j1] <- resj$resi
  
      bic2j[j1] <- resj$BIC2

      beta_arrayj2[,,j1] <- resj$betaest
      groupmatj2[,j1] <- resj$group
      etamatj2[,j1] <- resj$etaest
      groupetaj2[,j1] <- resj$groupc
    }
    
    
    indexmin <- which.min(bic2j)
    
    bic2k[j2] <- bic2j[indexmin]
    
    beta_arrayk2[,,j2] <- beta_arrayj2[,,indexmin]
    groupmatk2[,j2] <- groupmatj2[,indexmin]
    
    etamatk2[,j2] <- etamatj2[,indexmin]
    groupetak2[,j2]  <- groupetaj2[,indexmin]
  }
  

  indj <- which.min(bick2)
  bic2[k] <- bick2[indj]
  
  beta_array2[,,k] <- beta_arrayk2[,,indj]
  groupmat2[,k] <- groupmatk2[,indj]
  etamat2[,k] <- etamatk2[,indj]
  groupeta2[,k] <- groupetak2[,indj]
  
  print(j2)
}


res_c2 <-  refit_cohort2(year = year2, age = age2,y = dat$PropObese, x = x2,group.individual = groupmat2[,5],group.cohort = groupeta2[,5],model = "age")
