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


lam1 <- seq(0.05,0.2,by = 0.0025)
lam2 <- seq(0.001,0.004,by = 0.0005)
alp <- c(0.25,0.5,0.75,1,1.25,1.5)

bic2 <- rep(0, length(lam2))
beta_array2 <- array(0, dim = c(nage2,3,length(lam2)))
groupmat2 <- matrix(0, nage2, length(lam2))
etamat2 <- groupeta2 <- matrix(0, ncoh, length(lam2))


for(j2 in 1:length(lam2))
{
 
  bicj2 <- rep(0, length(alp))
  
  beta_arrayj2 <- array(0, dim = c(nage2,ncol(x2),length(alp)))
  groupmatj2 <- matrix(0, nage2, length(alp))
  
  etamatj2 <- matrix(0, ncoh, length(alp))
  groupetaj2  <- matrix(0, ncoh, length(alp))
  
  
  
  for(k in 1:length(alp))
  {
    betam2 <- betam02
    
    bic2k <- rep(0, length(lam1))
    
    beta_arrayk2 <- array(0, dim = c(nage2,ncol(x2),length(lam1)))
    groupmatk2 <- matrix(0, nage2, length(lam1))
    
    etamatk2 <- matrix(0, ncoh, length(lam1))
    groupetak2  <- matrix(0, ncoh, length(lam1))
    
    
    weightsk <- exp(alp[k]*(1-ordervec))
    
    for(j1 in 1:length(lam1))
    {
      resj <- Gr_cohort2(year = year2,age = age2,y = dat$PropObese,x = x2,
                         betam0 = betam2,model = "age",weights = weightsk,
                         lam1 = lam1[j1],lam2 = lam2[j2])
      
      
      betam2 <- resj$betaest
      
      bic2k[j1] <- resj$BIC2
      beta_arrayk2[,,j1] <- resj$betaest
      groupmat2k[,j1] <- resj$group
      etamatk2[,j1] <- resj$etaest
      groupetak2[,j1] <- resj$groupc
    }
    
    indexmin <- which.min(bic2k)
    bic2[k] <- bic2k[indexmin]
    bicc2[k] <- bicc2k[which.min(bicc2k)]
    beta_array2[,,k] <- beta_array2k[,,indexmin]
    groupmat2[,k] <- groupmat2k[,indexmin]
    groupmatc2[,k] <- groupmat2k[,which.min(bicc2k)]
  }
  
  
  
  indj <- which.min(bic2j)
  bic2[j2] <- bic2j[indj]
  
  beta_array1[,,j2] <- beta_arrayj[,,indj]
  groupmat1[,j2] <- groupmatj[,indj]
  etamat1[,j2] <- etamatj[,indj]
  groupeta1[,j2] <- groupetaj[,j2]
  
  print(j2)
}
