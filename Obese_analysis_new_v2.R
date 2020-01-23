library(Spgr)
library(plyr)
library(dplyr)
library(ggplot2)
library(mvtnorm)
#setwd("/Users/wangx172/Research/Obesity/")
#setwd("/Users/wangx172/Dropbox/Tanja&XinW/Rfiles/") 
source("Gr_cohort/refit_cohort.R")
source("Gr_cohort/Gr_cohort3.R")

dat <- read.csv("/Users/wangx172/Dropbox/Tanja&XinW/Newdata/AggObese1990_2017.csv")



dat2 <- arrange(dat, AGE, IYEAR)
dat2 <- filter(dat2, AGE !=80)
year2 <- dat2$IYEAR
age2 <- dat2$AGE
nage2 <- length(unique(age2))
y2 <- dat2$PropObese 
x2 <- cbind(1, scale(dat2$IYEAR),scale((dat2$IYEAR - mean(dat2$IYEAR))^2))
ncoh <- nage2 + length(unique(year2)) - 1
ages = unique(dat2$AGE)

########### first step for cohort ###
lamc <- seq(0.001,0.03,by = 0.0005)
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

groupc <- res_s1$groupc

######## step 2 #######


betam02 <- cal_initialrx(indexy = age2,y = scale(y2),x = x2)


lam2c <- seq(0.05,5,by = 0.05)
wmat <- matrix(0, nage2, nage2)
for(i in 1:(nage2-1))
  for(j in (i+1):(nage2))
    wmat[i,j] <- abs(i - j)
wmat <- wmat + t(wmat)
ordervec <- wmat[lower.tri(wmat)]

alp <- c(0, 0.25, 0.5, 0.75, 1, 1.25)

bic_s21 <- bic_s22 <-matrix(0, length(lam2c),length(alp))
group_s2 <- array(0, dim = c(nage2, length(lam2c),length(alp)))

betals = list()

for(k in 1:length(alp))
{
  weightsk <- exp(alp[k]*(1-ordervec))
  
  beta_arrayk = array(0, dim = c(nage2, 3, length(lam2c)))
  
  for(j in 1:length(lam2c))
  {
    res_s2j <- Gr_coef_only(year = dat2$IYEAR,age = age2, y = scale(y2),x = x2,betam0 = betam02,
                            ws = weightsk,model = "age", 
                            group.cohort = groupc,lam = lam2c[j],maxiter = 2000)
    bic_s21[j,k] <-res_s2j$BIC2
    bic_s22[j,k] <- res_s2j$BICc2
    group_s2[,j,k] <- res_s2j$group
    beta_arrayk[,,j] = res_s2j$betaest
  }
  
  betals[[k]] = beta_arrayk
}

inds21 = which(bic_s21==min(bic_s21),arr.ind = TRUE)

group_coef = group_s2[,10,4]


cbind(ages, group_coef)


betals[[4]][,,10]



plot(lam2c, betals[[4]][1,1,], type ="l",ylim = c(-2,1))
for(i in 2:62)
{
  lines(lam2c, betals[[4]][i,1,])
}


#### trace plot?? ######

ngroup4 = apply(group_s2[,,4],2,function(x){length(unique(x))})
plot(lam2c,ngroup4)
plot(bic_s21[,4],ngroup4)
plot(ngroup4[ngroup4 >=4 & ngroup4<=7], bic_s21[ngroup4 >=4 & ngroup4<=7,4])

bic_s21[ngroup4==6,4]
bic_s21[ngroup4==7,4]
bic_s21[ngroup4==5,4]
bic_s21[ngroup4==4,4]

c(min(bic_s21[ngroup4==4,4]),min(bic_s21[ngroup4==5,4]),min(bic_s21[ngroup4==6,4]),min(bic_s21[ngroup4==7,4]))

lam2c4 = lam2c[bic_s21[,4] == min(bic_s21[ngroup4==4,4])]
lam2c5 = lam2c[bic_s21[,4] == min(bic_s21[ngroup4==5,4])]
lam2c6 = lam2c[bic_s21[,4] == min(bic_s21[ngroup4==6,4])]
lam2c7 = lam2c[bic_s21[,4] == min(bic_s21[ngroup4==7,4])]
c(lam2c4, lam2c5, lam2c6, lam2c7)

### coefficients
betaest4 = unique(betals[[4]][,,lam2c == lam2c4])
betaest5 = unique(betals[[4]][,,lam2c == lam2c5])
betaest6 = unique(betals[[4]][,,lam2c == lam2c6])
betaest7 = unique(betals[[4]][,,lam2c == lam2c7])



summary(bic_s21[,4])
group_s2[,ngroup4==6,4]

group_coef4 = group_s2[,ngroup4==4,4][,which.min(bic_s21[ngroup4==4,4])]
group_coef5 = group_s2[,ngroup4==5,4][,which.min(bic_s21[ngroup4==5,4])]
group_coef7 = group_s2[,ngroup4==7,4][,which.min(bic_s21[ngroup4==7,4])]

#### refit 
uyear <- unique(dat2$IYEAR)
xm <- dat2$IYEAR - mean(uyear)
sigmax <- sum((uyear - mean(uyear))^2)/length(uyear)

x0 <- cbind(1, xm, xm^2 - sigmax)

res_fit <- refit_cohort2(year = year2,age = age2,y = dat2$PropObese,x = x0,group.individual = group_coef,group.cohort = groupc,model = "age")

res_fit2 <- refit_cohort2(year = year2,age = age2,y = scale(y2),x = x2,group.individual = group_coef,group.cohort = groupc,model = "age")

res_fitg4 <- refit_cohort2(year = year2,age = age2,y = scale(y2),x = x2,group.individual = group_coef4,group.cohort = groupc,model = "age")


save(betaest4, betaest5, betaest6, betaest7, res_fitg4, group_coef,
     group_coef4, group_coef5, group_coef7, file = "betaest4567.RData")

dfc <- data.frame(cohort = sort(unique((factor(year2 - age2)))), gc = groupc, ID = as.factor(groupc))
## cohort figure
pdf("doc/cohort_groupnew.pdf",height = 6,width = 8) # should be changed based on your local drive
ggplot(dfc,aes(x=cohort, y=gc, group=gc, shape= ID)) + geom_point() + 
  scale_shape_manual(values=1:nlevels(dfc$ID)) +
  labs(title = "", x="cohort", y="subgroups") +
  geom_point(size=2)+ theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_discrete("cohort", labels = seq(1911,1999,by = 5), breaks = seq(1911,1999,by = 5))+
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) + 
  scale_y_continuous("subgroups", labels = 1:5, breaks = 1:5)
dev.off()

## age figure 
ages = unique(dat2$AGE)
dfa <- data.frame(age = unique(age2), ga = group_coef, ID = as.factor(group_coef))
pdf("doc/age_groupnew.pdf",height = 6,width = 8) # should be changed based on your local drive
ggplot(dfa,aes(x=age, y=ga, group=ga, shape= ID)) + geom_point() + 
  scale_shape_manual(values=1:nlevels(dfa$ID)) +
  labs(title = "", x="ages", y="subgroups") +
  geom_point(size=2)+ theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("ages", labels = as.character(ages), breaks = ages)+
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous("subgroups", labels = as.character(1:6), breaks = 1:6)
dev.off()



###### age group with different number of groups #####
dfa2 <- data.frame(age = rep(unique(age2),4), ga = c(group_coef4, group_coef5, group_coef, group_coef7),
                   ID = as.factor(c(group_coef4, group_coef5, group_coef, group_coef7)), 
                   ngroup = rep(4:7, each = nage2)) 
dfa2 = dfa2 %>% mutate(
  ngroup = factor(ngroup, 4:7, labels = c("Model 1 (4 groups)", "Model 2 (5 groups)","Model 3 (6 groups)","Model 4 (7 groups)"))
)

pdf("doc/agegroups4fig.pdf",width = 6,height = 8)
ggplot(dfa2,aes(x=age, y=ga, group=ngroup, shape= ID)) + geom_point() + 
  scale_shape_manual("group",values=1:nlevels(dfa2$ID)) +
  labs(title = "", x="ages", y="subgroups") +
  geom_point(size=2)+ theme_bw() + 
  facet_wrap(~ngroup, nrow = 4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous("ages", labels = as.character(ages[seq(1,62,by=4)]), breaks = ages[seq(1,62,by=4)])+
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous("subgroups", labels = as.character(1:7), breaks = 1:7)
dev.off()
  
  
##### fitted curves ####
preddat <- dat2
preddat$est <- res_fit$estimates
preddat$group <- as.factor(rep(group_coef, each = length(unique(dat2$IYEAR))))
preddat$curve <- res_fit$curve

pdf("doc/group_curvesnew.pdf",height = 6,width = 8) 
ggplot(data = preddat) + 
  geom_point(aes(x = IYEAR, y = PropObese), alpha = 0.5) + 
  geom_line( aes(x = IYEAR, y = curve, linetype = group), size = 0.8) +
  ylab("obesity prevalence")+
  theme_bw()  +
  theme(legend.key.width=unit(1,"cm")) + 
  scale_x_continuous("year", labels = seq(1990,2017,by = 5), breaks = seq(1990,2017,by = 5))
dev.off()


save(bic_s11, groupc, bic_s21, group_s2, res_fit, file = "result_newv2.RData")




