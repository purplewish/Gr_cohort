#install.packages("devtools")
library(devtools)
#install_github("purplewish/Spgr")
setwd("~/Dropbox/Tanja&XinW")
#install.packages("plyr")
library(dplyr)
library(plyr)
#install.packages("Spgr")
#install.packages("Rcpp")
library(Spgr)
library(Rcpp)
obesedat <- read.csv("data/FullObeseYear1.csv")
xmean <- mean(unique(obesedat$AGE))
xvar <- var(unique(obesedat$AGE))
xmat <- cbind(1, scale(obesedat$AGE - xmean), scale((obesedat$AGE - xmean)^2 - xvar))
xmat1 <- cbind(1, scale(obesedat1$AGE - xmean), scale((obesedat1$AGE - xmean)^2 - xvar))
nyear <- length(unique(obesedat$IYEAR))

ddply(obesedat, .(IYEAR), summarize, nage = length(AGE))
?ddply

obesedat1 <- arrange(obesedat,desc(IYEAR))

betam01 <- cal_initialrx(indexy = obesedat$IYEAR,y = obesedat$PropObese,x = xmat,lam0 = 0.0001)
betam011 <- cal_initialrx(indexy = obesedat1$IYEAR,y = obesedat1$PropObese,x = xmat1,lam0 = 0.0001)

res1 <- Spgrrx(indexy = obesedat$IYEAR,y = obesedat$PropObese,x = xmat,weights = rep(1,nyear*(nyear-1)),betam0 = betam01,lam = 0.005)

res2 <- Spgrrx(indexy = obesedat1$IYEAR,y = obesedat1$PropObese,x = xmat1,weights = rep(1,nyear*(nyear-1)),betam0 = betam01,lam = 0.005)

help(Spgrrx)

lam1 <- seq(0.001,0.04,by = 0.0005)
bic1 <- rep(0, length(lam1))
beta_array <- array(0, dim = c(nyear,3,length(lam1)))
betam1 <- betam01
groupmat <- matrix(0, nyear, length(lam1))
ngroup <- rep(0, length(lam1))
for(j in 1:length(lam1))
{
  resj <- Spgrrx(indexy = obesedat$IYEAR,y = obesedat$PropObese,x = xmat,weights = rep(1,nyear*(nyear-1)),betam0 = betam01,lam = lam1[j])
  betam1 <- resj$beta
  bic1[j] <- BICcrx(obj = resj,indexy = obesedat$IYEAR,y = obesedat$PropObese,x = xmat,c0 = 1)
  beta_array[,,j] <- resj$beta
  ngroup[j] <- length(unique(resj$group))
  groupmat[,j] <- resj$group
}

plot(beta_array[,1,26])
plot(beta_array[,1,10])

groupnew <- data.frame(year = resj$index, group=as.factor(groupmat[,10]))
dat1 <- merge(obesedat,groupnew, by.x = "IYEAR",by.y = "year")
#install.packages("ggplot2")
#install.packages("cvTools")

library(ggplot2)
ggplot(dat1) + geom_point(aes(x = AGE, y = PropObese, color = group))


### cross validation based on nested ####
# 58 ages 
# 5 fold 
index <- as.numeric(cut(1:58, quantile(0:58,(0:5)/5)))
mse1 <- rep(0, length(lam1))
unage <- sort(unique(obesedat$AGE))
for(j in 1:length(lam1))
{
  msej <- 0
  for(j1 in 1:4)
  {
    index_train <- obesedat$AGE %in% unage[index<= j1]
    index_test <- obesedat$AGE %in% unage[index==j1+1]
    xmat_test <- xmat[index_test,]
    ytest <- obesedat$PropObese[index_test]
    
    resj1 <- Spgrrx(indexy = obesedat$IYEAR[index_train],y = obesedat$PropObese[index_train],x = xmat[index_train,],weights = rep(1,nyear*(nyear-1)),betam0 = betam01,lam = lam1[j])
    yearvj1 <- resj1$index
    year_test <- obesedat$IYEAR[index_test]
    betaj1 <- resj1$beta
    
    yhatj1 <- rep(0, length(ytest))
    
    for(j2 in 1:length(yearvj1))
    {
      indexj2 <- year_test== yearvj1[j2]
      yhatj1[indexj2] <- xmat_test[indexj2,] %*% betaj1[j2,]
      
    }
    
    msej <- msej + sum((yhatj1 - ytest)^2)
  }

  mse1[j] <- msej

}

mse1 <- mse1/58

### cross validation based on traditional cv###
yearv <- unique(obesedat$IYEAR)
library(cvTools)
set.seed(1236)
cvls <- list()
for(j2 in 1:nyear)
{
  temp <- cvFolds(length(unage),K = 5)
  cvls[[j2]] <- cbind(yearv[j2],temp$which, temp$subsets)
}


mse11 <- rep(0, length(lam1))
for(j in 1:length(lam1))
{
  msej <- 0
  for(j1 in 1:5)
  {
    index_train <- rep(0,nrow(obesedat))
    for(j2 in 1:nyear)
    {
      index_train[obesedat$IYEAR==yearv[j2]][cvls[[j2]][cvls[[j2]][,2]!=j1,3]] <- 1
    }
    
    index_test <- 1 - index_train
    
    index_test <- index_test==1
    index_train <- index_train==1
    
    xmat_test <- xmat[index_test,]
    ytest <- obesedat$PropObese[index_test]
    
    resj1 <- Spgrrx(indexy = obesedat$IYEAR[index_train],y = obesedat$PropObese[index_train],x = xmat[index_train,],weights = rep(1,nyear*(nyear-1)),betam0 = betam01,lam = lam1[j])
    yearvj1 <- resj1$index
    year_test <- obesedat$IYEAR[index_test]
    betaj1 <- resj1$beta
    
    yhatj1 <- rep(0, length(ytest))
    
    for(j2 in 1:length(yearvj1))
    {
      indexj2 <- year_test== yearvj1[j2]
      yhatj1[indexj2] <- xmat_test[indexj2,] %*% betaj1[j2,]
      
    }
    
    msej <- msej + sum((yhatj1 - ytest)^2)
  }
  
  mse11[j] <- msej
}

res11 <- Spgrrx(indexy = obesedat$IYEAR,y = obesedat$PropObese,x = xmat,weights = rep(1,nyear*(nyear-1)),betam0 = betam01,lam = 0.003)


