#This is for exercise 1, Bootstrapping section B, problem 1

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


#B.1) simulate multivariate normal random variables given mu/sigma 
mu = c(5,0)
s = matrix(c(2,1,1,2), 2, 2, byrow = T)
N = 1000
sim.nums = draw.multiVarNorm(mu,s)
#plot(sim.nums[,1],sim.nums[,2],xlim=c(-5,5),ylim=c(-5,5))

