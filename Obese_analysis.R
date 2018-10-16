library(Spgr)
library(plyr)
library(ggplot2)

dat <- read.csv("Dropbox/Tanja&XinW/Rfiles/CBD-O/FullObeseYear1.csv")
dat <- arrange(dat, IYEAR, AGE)
year <- dat$IYEAR
nyear <- length(unique(year))
age <- dat$AGE
y <- dat$PropObese
x <- cbind(1, scale(dat$AGE), scale(dat$AGE^2))
nu <- 1
gam <- 3
weights <- rep(1, nyear*(nyear-1))
betam01 <- cal_initialrx(indexy = dat$IYEAR,y = y,x = x) # from Spgr package 


###### find group for years #####
lam1 <- seq(0.001,0.04,by = 0.0005)
bic1 <- rep(0, length(lam1))
beta_array <- array(0, dim = c(nyear,3,length(lam1)))
betam1 <- betam01
groupmat <- matrix(0, nyear, length(lam1))
ngroup <- rep(0, length(lam1))
for(j in 1:length(lam1))
{
  resj <- Gr_cohort(year = year,age = age,y = y,x = x,betam0 = betam1,model = "year",weights = weights,lam = lam1[j])
  betam1 <- resj$betaest
  bic1[j] <- resj$BIC
  beta_array[,,j] <- resj$betaest
  ngroup[j] <- length(unique(resj$group))
  groupmat[,j] <- resj$group
}

plot(beta_array[,,29])

cbind(unique(year),groupmat[,29])
unique(beta_array[,,29])
res_year <- refit_cohort(year = year,age = age,x = x,group = groupmat[,29],model = "year") # refit with estimated group structure. 


preddat <- dat
preddat$est <- res_year$estimates
preddat$group <- as.factor(rep(groupmat[,29], each = length(unique(age))))

ggplot(data = preddat, aes(x = AGE, y = PropObese, group = group, color = group)) + 
  geom_point() + theme_bw()


#### find groups for ages  need to discuss #####
library(splines)
dat2 <- arrange(dat, AGE, year)
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
  beta_array2k <- array(0, dim = c(nage2,ncol(x2),length(lam2)))
  groupmat2k <- matrix(0, nage2, length(lam2))
  weightsk <- exp(alp[k]*(1-ordervec))
  
  for(j in 1:length(lam2))
  {
    resj <- Gr_cohort(year = year2,age = age2,y = y2,x = x2,betam0 = betam2,model = "age",
                      weights = weightsk,lam = lam2[j])
    betam2 <- resj$betaest
    bic2k[j] <- resj$BIC
    beta_array2k[,,j] <- resj$betaest
    groupmat2k[,j] <- resj$group
  }
  
  indexmin <- which.min(bic2k)
  bic2[k] <- bic2k[indexmin]
  beta_array2[,,k] <- beta_array2k[,,indexmin]
  groupmat2[,k] <- groupmat2k[,indexmin]
}


