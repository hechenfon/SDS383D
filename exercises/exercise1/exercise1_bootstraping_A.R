#This is for exercise 1, Bootstrapping section, A

boot.raw<-function(X,y,n = 1e4) {
  #bootstrapping on raw sampleto calculate the estimator variance
  len = dim(X)[1]
  p = dim(X)[2]
  beta.all<-matrix(nrow=n,ncol=p)
  for (i in c(1:n)) {
    id.boot = sample(c(1:len),len,replace=TRUE)
    X.new = X[id.boot,]
    y.new = y[id.boot]
    betahat = solve(crossprod(X.new)) %*% t(X.new) %*% y.new
    beta.all[i,] = betahat
  }

  boot.res = cov(beta.all)
  return (boot.res)
}

# A). main starts, estimate sigma variance
library(mlbench)
ozone = data(Ozone, package='mlbench')
ozone = na.omit(Ozone)[,4:13]
y = ozone[,1]
x = as.matrix(ozone[,2:10])
x = cbind(1,x)
#estimate variance based on bootstrapping on raw samples
betacov.res = boot.raw(x,y)  
