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

sig2est = 0.083

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

for(j in 1:length(lam2c))
{
  res_s4j = Spgrrx(indexy = age2, y = ysim,x = x2,weights = weights4,
                           betam0 = betam04,lam = lam2c[j], maxiter = 1000)
  bic_s4[j] <-  ###calculate BIC
}



#### refit ##
res_rfit4 = refit_group(year = dat2$IYEAR, age = age2, y = ysim, x = x2, group = group_coef4, model = "age")
betaestr = matrix(res_rfit4$betaest, ncol = 3, byrow = TRUE)

norm(betaest40 - betaestr,"f")
norm(betaest40 - res_s4j$betaest,"f")



