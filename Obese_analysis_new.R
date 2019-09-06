library(Spgr)
library(plyr)
library(dplyr)
library(ggplot2)
#setwd("/Users/wangx172/Research/Obesity/")
setwd("/Users/wangx172/Dropbox/Tanja&XinW/Rfiles/") 
source("Gr_cohort/refit_cohort.R")
source("Gr_cohort/Gr_cohort3.R")

dat <- read.csv("/Users/wangx172/Dropbox/Tanja&XinW/Newdata/AggObese1990_2017.csv")


########### first step for cohort ###
dat2 <- arrange(dat, AGE, IYEAR)
dat2 <- filter(dat2, AGE !=80)
year2 <- dat2$IYEAR
age2 <- dat2$AGE
nage2 <- length(unique(age2))
y2 <- dat2$PropObese 
x2 <- cbind(1, scale(dat2$IYEAR),scale((dat2$IYEAR - mean(dat2$IYEAR))^2))
ncoh <- nage2 + length(unique(year2)) - 1



lamc <- seq(0.001,0.02,by = 0.00025)
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
                         model = "age", group.individual = 1:nage2,lam2 = lamc[which.min(bic_s11)],maxiter = 2000)

res_s1 <- Gr_cohort_only(year = dat2$IYEAR,age = age2, y = scale(y2),x = x2,
                         model = "age", group.individual = group_coef,lam2 = lamc[6],maxiter = 2000)


groupc <- res_s1$groupc


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


lam2c <- exp(seq(-2, 6, by = 0.05))
bic_s21 <- rep(0, length(lam2c))
beta_array_s2 <- array(0, dim = c(nage2,3,length(lam2c)))
group_s2 <- matrix(0, nage2, length(lam2c))

weightsk <- exp(alp[5]*(1-ordervec))
betam02j <- betam02

for(j in 1:length(lam2c))
{
  res_s2j <- Gr_coef_only(year = dat2$IYEAR,age = age2, y = scale(y2),x = x2,betam0 = betam02,
                          ws = weightsk,
                          model = "age", group.cohort = groupc,lam = lam2c[j],maxiter = 2000)
  betam02j <- res_s2j$betaest
  bic_s21[j] <-res_s2j$BIC2
  beta_array_s2[,,j] <- res_s2j$betaest
  group_s2[,j] <- res_s2j$group
}


res_s2 <- Gr_coef_only(year = dat2$IYEAR,age = age2, y = scale(y2),x = x2,betam0 = betam02,
                        ws = weightsk,
                        model = "age", group.cohort = groupc,lam = 100,maxiter = 2000)

#### refit 
uyear <- unique(dat2$IYEAR)
xm <- dat2$IYEAR - mean(uyear)
sigmax <- sum((uyear - mean(uyear))^2)/length(uyear)

x0 <- cbind(1, xm, xm^2 - sigmax)

res_fit <- refit_cohort2(year = year2,age = age2,y = dat2$PropObese,x = x0,group.individual = group_coef,group.cohort = groupc,model = "age")

save(group_coef, groupc, res_fit, file = "result_new.RData")

load("result_new.RData")

dfc <- data.frame(cohort = sort(unique((factor(year2 - age2)))), gc = groupc, ID = as.factor(groupc))
## cohort figure
pdf("doc/cohort_group.pdf",height = 6,width = 8) # should be changed based on your local drive
ggplot(dfc,aes(x=cohort, y=gc, group=gc, shape= ID)) + geom_point() + 
  scale_shape_manual(values=1:nlevels(cc$ID)) +
  labs(title = "", x="cohort", y="subgroups") +
  geom_point(size=2)+ theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_discrete("cohort", labels = seq(1911,1999,by = 5), breaks = seq(1911,1999,by = 5))+
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) + 
  scale_y_continuous("subgroups", labels = 1:7, breaks = 1:7)
dev.off()

## age figure 
dfa <- data.frame(age = unique(age2), ga = group_coef, ID = as.factor(group_coef))
## cohort figure
pdf("doc/age_group.pdf",height = 6,width = 8) # should be changed based on your local drive
ggplot(dfa,aes(x=age, y=ga, group=ga, shape= ID)) + geom_point() + 
  scale_shape_manual(values=1:nlevels(cc$ID)) +
  labs(title = "", x="ages", y="subgroups") +
  geom_point(size=2)+ theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("ages", labels = as.character(ages), breaks = ages)+
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous("subgroups", labels = as.character(1:7), breaks = 1:7)
dev.off()


preddat <- dat2
preddat$est <- res_fit$estimates
preddat$group <- as.factor(rep(group_coef, each = length(unique(dat2$IYEAR))))
preddat$curve <- res_fit$curve

pdf("doc/group_curves.pdf",height = 6,width = 8) 
ggplot(data = preddat) + 
  geom_point(aes(x = IYEAR, y = PropObese), alpha = 0.5) + 
  geom_line( aes(x = IYEAR, y = curve, linetype = group), size = 0.8) +
  ylab("obesity prevalence")+
  theme_bw()  +
  theme(legend.key.width=unit(1,"cm")) + 
  scale_x_continuous("year", labels = seq(1990,2017,by = 5), breaks = seq(1990,2017,by = 5))
dev.off()

ggplot(data = preddat) + 
  geom_point(aes(x = IYEAR, y = PropObese), alpha = 0.5) + 
  geom_line( aes(x = IYEAR, y = curve, linetype = group)) +
  theme_bw()  +
  theme()
  facet_wrap(~group)
  scale_x_continuous("year", labels = seq(1990,2017,by = 5), breaks = seq(1990,2017,by = 5))

####### diagnois plot #####
library(ggplot2)
library(dplyr)
dat3 = dat2
dat3$group = findInterval(dat3$AGE,vec = c(0,30,40,50,60,70,80),rightmost.closed = TRUE)
ggplot(data = filter(dat3,AGE >=23),aes(x = IYEAR, y = PropObese, group = AGE, color = AGE)) + geom_line() + facet_wrap(~group)


# dat0 = merge(dat, dat2, by = c("IYEAR","AGE"),all = TRUE)
# 
# dat4 = dat
# dat4$group = findInterval(dat4$AGE,vec = c(0,30,40,50,60,70,80),rightmost.closed = TRUE)
# ggplot(data = dat4,aes(x = IYEAR, y = PropObese, group = AGE, color = AGE)) + geom_line() + facet_wrap(~group)

plot(dat3$PropObese[dat3$AGE==80]
,dat4$PropObese[dat4$AGE==80])
abline(0,1)
