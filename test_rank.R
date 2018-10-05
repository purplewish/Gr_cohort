ny <- 4 #number of ages considered
na <- 8 #number of years considered
n <- na * ny 
nc <- na + ny - 1

group <- c(1,1,2,3)

#####X1   #this is Xa in Currie(2016)
X1 <- kronecker(diag(ny), rep(1, na)) 
X1
dim(X1)

#####X2  

x <- 1:(1+na-1) #ages
xm <- x 
X2 <- kronecker(diag(ny), xm)

######X3
sigmax <- sum(xm^2)/na
xq <- xm^2 
X3 <- kronecker(diag(ny), xq)
dim(X3)
######Xc
cohort = c()
for(i in 1:ny){
  cohort <- cbind(cohort, (na:1+i-1))
}


Xc = c()
for(j in 1:ny){
  for(i in 1:na){
    Xc <- rbind(Xc, as.numeric((1:nc) == cohort[i, j]))
  }
}
dim(Xc)
################################M7#######################
X = cbind(X1, X2, X3, Xc)
y1 = rowSums(X3)
x1 <- cbind(X2, Xc)
round(Solve(x1,y1),2)

qr(X)$rank

X_v1 <- X
X_v1[,1:12] <- X_v1[,c(1,5,9,2,6,10,3,7,11,4,8,12)]

ng <- length(unique(group))
W <- matrix(0, ny, ng)
W[cbind(1:ny,group)] <- 1
W <- W %x% diag(1,3)
Ux <- cbind(X_v1[,1:(ncol(X1)*3)]%*%W,Xc)

qr(Ux)$rank
