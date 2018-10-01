############# refit algorithm #######
# refit based on a given group information

# only 
refit_cohort <- function(group,nr,x, y, indexy)
{
  ncx <- ncol(x)
  Xm <- matrix(0, nrow(x), nr*ncx)
  ym <- rep(0, nrow(x))
  uniqxy <- unique(indexy)
  nJ <- rep(0, nr)
  for(i in 1:nr)
  {
    nJ[i] <- sum(indexy == uniqxy[i])
    Xm[indexy == uniqxy[i],(ncx*(i-1) + 1) : (ncx*i)] <- x[indexy == uniqxy[i],]/sqrt(nJ[i])
    ym[indexy == uniqxy[i]] <- y[indexy == uniqxy[i]]/sqrt(nJ[i])
  }
  
  ng <- length(unique(group))
  W <- matrix(0, nr, ng)
  W[cbind(1:nr,group)] <- 1
  W <- W %x% diag(1,ncx)
  Ux <- Xm%*%W
  est <- solve(t(Ux)%*%Ux)%*%t(Ux)%*%ym
  
  sig2 <- sum((ym - Ux%*%est)^2)/nr
  
  out <- list(est = est, sig2 = sig2)
  return(out)
}


