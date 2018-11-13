library(Spgr)
library(plyr)
library(ggplot2)
setwd("Research/Obesity/")
source("Gr_cohort/Gr_cohort.R")
source("Gr_cohort/refit_cohort.R")
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

y <- dat$PropObese
x <- cbind(1, scale(dat$AGE), scale(dat$AGE^2))
nu <- 1
gam <- 3
weights <- rep(1, nyear*(nyear-1))
betam01 <- cal_initialrx(indexy = dat$IYEAR,y = y,x = x) # from Spgr package 


lam1 <- seq(0.001,0.04,by = 0.0005)
bic1 <- rep(0, length(lam1))
beta_array <- array(0, dim = c(nyear,3,length(lam1)))
betam1 <- betam01
groupmat <- matrix(0, nyear, length(lam1))
ngroup <- rep(0, length(lam1))
resi1 <- rep(0, length(lam1))
for(j in 1:length(lam1))
{
  resj <- Gr_cohort(year = year,age = age,y = y,x = x,betam0 = betam1,model = "year",weights = weights,lam = lam1[j])
  betam1 <- resj$betaest
  bic1[j] <- resj$BIC
  resi1[j] <- resj$resi
  beta_array[,,j] <- resj$betaest
  ngroup[j] <- length(unique(resj$group))
  groupmat[,j] <- resj$group
}

which.min(bic1)
plot(beta_array[,,29])

cbind(unique(year),groupmat[,29])
unique(beta_array[,,29])


### original scale 
xm <- dat$AGE - mean(uage)
sigmax <- sum((uage - mean(uage))^2)/length(uage)

x0 <- cbind(1, xm, xm^2 - sigmax)
res_year <- refit_cohort(year = year,age = age, y = dat$PropObese,  x = x0,group = groupmat[,29],model = "year") # refit with estimated group structure. 

res_year2 <- refit_cohort2(year = year,age = age, y = dat$PropObese,x = x0,group.individual = groupmat[,29], group.cohort = c(rep(1:42,each = 2),43),model = "year") 

res_gr_year <- refit_group(year = year, age = age, y = dat$PropObese, x = x0,group = groupmat[,29],model="year") ## fit linear model without cohort 



preddat <- dat
preddat$est <- res_year$estimates
preddat$group <- as.factor(rep(groupmat[,29], each = length(unique(age))))
preddat$curve <- res_year$curve
preddat$curve_gr <- res_gr_year$estimates

ggplot(data = preddat) + 
  geom_point(aes(x = AGE, y = PropObese, group = group, color = group)) + 
  geom_line(aes(x = AGE, y = curve, group = group, color = group)) +
 # geom_line(aes(x = AGE, y = curve_gr, group = group, color = group)) +
  theme_bw() 





####### groups for age #######
### based on quadratic #######

dat2 <- arrange(dat, AGE, year)
year2 <- dat2$IYEAR
age2 <- dat2$AGE
nage2 <- length(unique(age2))
y2 <- dat2$PropObese 
x2 <- cbind(1, scale(dat2$IYEAR),scale(dat2$IYEAR^2))


# weights for age is defined based on the age difference 
wmat <- matrix(0, nage2, nage2)
for(i in 1:(nage2-1))
  for(j in (i+1):(nage2))
    wmat[i,j] <- abs(i - j)
wmat <- wmat + t(wmat)
ordervec <- wmat[lower.tri(wmat)]

betam02 <- cal_initialrx(indexy = age2,y = y2,x = x2) # from Spgr package 


lam2 <- seq(0.05,0.2,by = 0.0025)
alp <- c(0.25,0.5,0.75,1,1.25,1.5)
bic2 <- rep(0, length(alp))
bicc2 <- rep(0, length(alp))
beta_array2 <- array(0, dim = c(nage2,ncol(x2),length(alp)))
groupmat2 <- matrix(0, nage2, length(alp))
groupmatc2 <- matrix(0, nage2, length(alp))
ngroup2 <- rep(0, length(alp))

for(k in 1:length(alp))
{
  betam2 <- betam02
  bic2k <- rep(0, length(lam2))
  bicc2k <- rep(0, length(lam2))
  resik <- rep(0, length(lam2))
  beta_array2k <- array(0, dim = c(nage2,ncol(x2),length(lam2)))
  groupmat2k <- matrix(0, nage2, length(lam2))
  weightsk <- exp(alp[k]*(1-ordervec))
  
  for(j in 1:length(lam2))
  {
    resj <- Gr_cohort(year = year2,age = age2,y = y2,x = x2,betam0 = betam2,model = "age",
                      weights = weightsk,lam = lam2[j])
    betam2 <- resj$betaest
    bic2k[j] <- resj$BIC
    bicc2k[j] <- resj$BICc
    resik[j] <- resj$resi
    beta_array2k[,,j] <- resj$betaest
    groupmat2k[,j] <- resj$group
  }
  
  indexmin <- which.min(bic2k)
  bic2[k] <- bic2k[indexmin]
  bicc2[k] <- bicc2k[which.min(bicc2k)]
  beta_array2[,,k] <- beta_array2k[,,indexmin]
  groupmat2[,k] <- groupmat2k[,indexmin]
  groupmatc2[,k] <- groupmat2k[,which.min(bicc2k)]
}


indmin2 <- which.min(bic2)

cbind(unique(age2),groupmat2[,indmin2])

res_age<- refit_cohort(year = year2,age = age2,y = dat2$PropObese,  x = x2,group = groupmat2[,indmin2],model = "age") # refit with estimated group structure. 

preddat2 <- dat2
preddat2$est <- res_age$estimates
preddat2$group <- as.factor(rep(groupmat2[,indmin2], each = length(unique(year2))))
preddat2$curve <- res_age$curve

ggplot(data = preddat2) + 
  geom_point(aes(x = IYEAR, y = PropObese, group = group, color = group)) + 
  geom_line(aes(x = IYEAR, y = curve, group = group, color = group))+
  theme_bw() 



#### nonparametric  find groups for ages  need to discuss  #####

library(splines)
year2 <- dat2$IYEAR
age2 <- dat2$AGE
nage2 <- length(unique(age2))
y2 <- dat2$PropObese

#x2 <- cbind(1, scale(dat2$IYEAR), scale(dat2$IYEAR^2))
#x2 <- cbind(1, scale(dat2$IYEAR))
x2 <- bs(x = dat2$IYEAR, intercept = TRUE, knots = quantile(unique(year2),c(0.25,0.5,0.75)),Boundary.knots = range(unique(year2)))
nu <- 1
gam <- 3

# weights for age is defined based on the age difference 
wmat <- matrix(0, nage2, nage2)
for(i in 1:(nage2-1))
  for(j in (i+1):(nage2))
    wmat[i,j] <- abs(i - j)
wmat <- wmat + t(wmat)

ordervec <- wmat[lower.tri(wmat)]

betam02 <- cal_initialrx(indexy = age2,y = y2,x = x2) # from Spgr package 


lam2 <- seq(0.005,0.06,by = 0.0025)
alp <- c(0.1,0.3,0.5,0.7,1)
bic2 <- rep(0, length(alp))
beta_array2 <- array(0, dim = c(nage2,ncol(x2),length(alp)))
groupmat2 <- matrix(0, nage2, length(alp))
ngroup2 <- rep(0, length(alp))

for(k in 1:length(alp))
{
  betam2 <- betam02
  bic2k <- rep(0, length(lam2))
  resik <- rep(0, length(lam2))
  beta_array2k <- array(0, dim = c(nage2,ncol(x2),length(lam2)))
  groupmat2k <- matrix(0, nage2, length(lam2))
  weightsk <- exp(alp[k]*(1-ordervec))
  
  for(j in 1:length(lam2))
  {
    resj <- Gr_cohort(year = year2,age = age2,y = y2,x = x2,betam0 = betam2,model = "age",
                      weights = weightsk,lam = lam2[j])
    betam2 <- resj$betaest
    bic2k[j] <- resj$BIC
    resik[j] <- resj$resi
    beta_array2k[,,j] <- resj$betaest
    groupmat2k[,j] <- resj$group
  }
  
  indexmin <- which.min(bic2k)
  bic2[k] <- bic2k[indexmin]
  beta_array2[,,k] <- beta_array2k[,,indexmin]
  groupmat2[,k] <- groupmat2k[,indexmin]
}


