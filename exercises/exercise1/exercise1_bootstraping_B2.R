# this is for exercise 1, bootstrapping section, problem B) 2)
# here to use the public available package 'mvnmle'

library('mvnmle')

#generate sample data to test
draw.multiVarNorm<-function(mu,sigma,n=1e3) {
  #simulate random variables based on given mean vector and covariance matrix
  dec = eigen(sigma)
  D_half = sqrt(diag(dec$values))
  V = dec$vectors
  L = V %*% D_half
  pre.trans = matrix(rnorm(length(mu)*n,mean=0,sd=1),nrow = n, byrow=TRUE)
  aft.trans = L %*% t(pre.trans)+mu
  return(t(aft.trans))
}
mu = c(5,0)
s = matrix(c(2,1,1,2), 2, 2, byrow = T)
N = 1000
sim.nums = draw.multiVarNorm(mu,s)
#generate sample data ending

#MLE estimation of mu/sigma based on sample data
est = mlest(sim.nums)

mu.hat = est$muhat
sigma.hat = est$sigmahat


