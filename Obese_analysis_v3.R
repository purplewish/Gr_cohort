####### two steps ####
library(Spgr)
library(plyr)
library(dplyr)
library(ggplot2)
setwd("Research/Obesity/")
source("Gr_cohort/Gr_cohort.R")
source("Gr_cohort/refit_cohort.R")
source("Gr_cohort/Gr_cohort2.R")
dat <- read.csv("/Users/wangx172/Dropbox/Tanja&XinW/Rfiles/CBD-O/FullObeseYear1.csv")
source("Gr_cohort/Gr_cohort3.R")


########### first step for cohort ###
dat2 <- arrange(dat, AGE, IYEAR)
year2 <- dat2$IYEAR
age2 <- dat2$AGE
nage2 <- length(unique(age2))
y2 <- dat2$PropObese 
x2 <- cbind(1, scale(dat2$IYEAR),scale((dat2$IYEAR - mean(dat2$IYEAR))^2))
ncoh <- nage2 + length(unique(year2)) - 1

lamc <- seq(0.01,0.04,length = 30)
bic_s11 <- bic_s12 <- rep(0,length(lamc))

for(j in 1:length(lamc))
{
  resj <- Gr_cohort_only(year = dat2$IYEAR,age = age2, y = scale(y2),x = x2,
                         model = "age", group.individual = 1:nage2,lam2 = lamc[j],maxiter = 2000)
  bic_s11[j] <- resj$BIC2
  bic_s12[j] <- resj$BICc2
}

which.min(bic_s11)
res_s1 <- Gr_cohort_only(year = dat2$IYEAR,age = age2, y = scale(y2),x = x2,
                         model = "age", group.individual = 1:nage2,lam2 = lamc[9],maxiter = 2000)

groupc <- res_s1$groupc



### second step for regression coefficients ##
lam2c <- seq(0.02, 0.4, by = 0.02)
betam02 <- cal_initialrx(indexy = age2,y = scale(y2),x = x2)

wmat <- matrix(0, nage2, nage2)
for(i in 1:(nage2-1))
  for(j in (i+1):(nage2))
    wmat[i,j] <- abs(i - j)
wmat <- wmat + t(wmat)
ordervec <- wmat[lower.tri(wmat)]

alp <- c(0, 0.25, 0.5, 0.75, 1, 1.25)

bic_s21a <- bic_s22 <- rep(0, length(alp))
beta_array_s2a <- array(0, dim = c(nage2,3,length(alp)))
group_s2a <- matrix(0, nage2, length(alp))

for(k in 1:length(alp))
{
  bic_s21 <- bic_s22 <- rep(0, length(lam2c))
  beta_array_s2 <- array(0, dim = c(nage2,3,length(lam2c)))
  group_s2 <- matrix(0, nage2, length(lam2c))
  
  weightsk <- exp(alp[k]*(1-ordervec))
  betam02j <- betam02
  
  for(j in 1:length(lam2c))
  {
    res_s2j <- Gr_coef_only(year = dat2$IYEAR,age = age2, y = scale(y2),x = x2,betam0 = betam02j,
                            ws = weightsk,
                            model = "age", group.cohort = groupc,lam = lam2c[j],maxiter = 2000)
    betam02j <- res_s2j$betaest
    bic_s21[j] <-res_s2j$BIC2
    bic_s22[j] <- res_s2j$BICc2
    beta_array_s2[,,j] <- res_s2j$betaest
    group_s2[,j] <- res_s2j$group
  }
  
  inds2 <- which.min(bic_s21)
  beta_array_s2a[,,k] <- beta_array_s2[,,inds2]
  group_s2a[,k] <- group_s2[,inds2]
  bic_s21a[k] <- min(bic_s21)
  
}

ind2 <- which.min(bic_s21a)
group_coef <- group_s2a[,ind2]


uyear <- unique(dat2$IYEAR)
xm <- dat2$IYEAR - mean(uyear)
sigmax <- sum((uyear - mean(uyear))^2)/length(uyear)

x0 <- cbind(1, xm, xm^2 - sigmax)

res_fit <- refit_cohort2(year = year2,age = age2,y = dat2$PropObese,x = x0,group.individual = group_coef,group.cohort = groupc,model = "age")

save(group_coef, groupc, res_fit, file = "result_v2.RData")


preddat <- dat2
preddat$est <- res_fit$estimates
preddat$group <- as.factor(rep(group_coef, each = length(unique(dat2$IYEAR))))
preddat$curve <- res_fit$curve


ggplot(data = preddat) + 
  geom_point(aes(x = IYEAR, y = PropObese), alpha = 0.5) + 
  geom_line( aes(x = IYEAR, y = curve, color = group)) +
  #geom_line(aes(x = IYEAR, y = est, group = group, color = group)) +
  theme_bw() 

age_df <- data.frame(age = unique(age2),group = group_coef)
plot(age_df)
