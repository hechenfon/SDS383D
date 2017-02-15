
library(ggplot2)
library(mvtnorm)

# sample beta
rbeta = function(m2,k2,w2){
  rbeta = sample_mvn(t(m2), solve(k2)/w2, 1)
  return(rbeta)
}

rlambda<-function(h, n, beta, w, Y, X) {
  #sampling from conditional distribution of lambda given y, beta and w
  #a = (h+1)/2
  #b = (1/2)*(w*rss+h)
  #lambda = rgamma(n, shape = a, rate = b)
  #lambda = rgamma(n, shape = rep(a,n), rate = b)
  lambda = rgamma(n, rep(h/2+1/2,n),h/2+w*(Y-X %*% beta)^2/2)
  #return(lambda)
  
  
} 

rw<-function(d2,eta2) {
  #sampling from conditional distribution of w given y, lambda
  w = rgamma(1, shape = d2/2, rate = eta2/2)
  return(w)
}

sample_mvn<-function(mu,cov,n=1e3) {
  #simulate random variables based on given mean vector and covariance matrix
  p = length(mu)
  unit_norm_sample = matrix(nrow = n, ncol = p)
  mvn_sample = matrix(nrow = n, ncol = p)
  for (j in 1:p){
    unit_norm_sample[,j] = rnorm(n)
  }
  L = chol(cov)
  mvn_sample =  t(L%*% t(unit_norm_sample)) + mu
  return(mvn_sample)
}



gdp<-read.csv("gdpgrowth.csv",header=T,stringsAsFactors = FALSE)

X = cbind(rep(1,dim(gdp)[1]),gdp$DEF60) #add one column for interception
Y = gdp$GR6096

#######BAYESIAN LINEAR REGRESSION#################
#define the prior parameters
k = matrix(data=c(0.001,0,0,0.001),nrow = 2, ncol = 2,byrow = T)  #diagnal & pretty vague precision matrix for w
d = 1
eta = 0.001
m = rep(0,dim(k)[1])
#calculate posterior distribution for w & beta
n = dim(X)[1]
p = dim(X)[2]
#d*
d2 = n+d
#k*
k2 = t(X) %*% X+k
#m*
m2 = solve(k2) %*% (t(X) %*% Y+k %*% m)
#eta*
eta2 = eta+t(Y) %*% Y+t(m) %*% k %*% m-t(m2) %*% k2 %*% m2

#posterior distribution of beta: p(beta|y)
#according to (C), p(beta|y) is a t-distribution: t(m2,d2,(eta2/d2)K2), posterior estimate of beta is m2
beta_post= m2

#pred_y = X %*% beta_post










##############HEAVY-TAILED part###############
iter = 1e4
x_dim = 2
p_dim = 2
d = 0.01
eta = 0.01
k = matrix(data=c(0.001,0,0,0.001),nrow = 2, ncol = 2,byrow = T) 
m = rep(0,dim(k)[1])
lambda = c(0, 0)
n = dim(X)[1]
h = 0.0001

beta = matrix(nrow = iter, ncol = p_dim)
w = rep(0,n)
#w[1] = rw(d, eta)
w[1] = 1
lambda = matrix(nrow = iter, ncol = n)
# s
d2 = n+d
#k*
lambda_p = diag(rep(1,n))

k2 = t(X) %*%  lambda_p %*% X+k
#m*
m2 = solve(k2) %*% (t(X) %*% lambda_p %*% Y+k %*% m)
#eta*
eta2 = eta+t(Y) %*% lambda_p %*% Y+t(m) %*% k %*% m-t(m2) %*% k2 %*% m2


## gibs sampler##
for(i in 2:iter){
  #beta[i,] = rbeta(m2, k2, w[i-1])
  beta[i,] = rmvnorm(1, m2, sigma = solve(k2)/w[i-1])
  
  # update beta
  
  #beta[i,] = rbeta(y, w, lambda, x, k, m)
  w[i] =rw(d2, eta2)
  
  #rss
  #rss = crossprod(Y - X%*%beta[i,])
  # sample each lambda 
  #lambda[i,] = rlambda(h, n, beta[i,], w[i], rss)
  lambda[i,] = rlambda(h, n, beta[i,], w[i], Y,X)
  #lambda[i,] = rgamma(n, rep(h/2+1/2,n),h/2+w[i]*(Y-X %*% beta[i,])^2/2)
  lambda_p = diag(lambda[i,])
  k2 = t(X) %*%  lambda_p %*% X+k
  m2 = solve(k2) %*% (t(X) %*% lambda_p %*% Y+k %*% m)
  eta2 = eta+t(Y) %*% lambda_p %*% Y+t(m) %*% k %*% m-t(m2) %*% k2 %*% m2

}

beta.avg = colMeans(beta[(iter-1000):iter,])

#############for plotting##########################
plot(X[,2],Y,xlab='DEF60',ylab='GR6096')
abline(a=beta.avg[1],b=beta.avg[2])
abline(a=beta_post[1],b=beta_post[2],col='red')

legend('bottomright',cex=.8,legend=c('heavy-tailed','not-heavy-tailed'),col=c('black','red'),lty=1)
