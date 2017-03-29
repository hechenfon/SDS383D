#Exercise 4, Problem 1
library('MCMCpack')
library(mvtnorm)

rm(list=ls())

genInd <- function(x, t) {
  M = matrix(0,length(x),length(t))
  for (i in c(1:length(x))) {
    yy = x[i]==t
    M[i,] = yy
  }
  return(M*1)
}


math<-read.csv('../../data/mathtest.csv')


#mean of each group
y.bar = aggregate(math,by=list(math$school),FUN=mean)
#number within each group
sample.num = as.vector(table(math$school))
#mean value of each group
plot(sample.num,y.bar$mathscore)  #high/low mean tend to be groups with small samples

#fit hierachical modelling, with prior:
#mu: prop to 1
#sigmaS: prop to 1/sigma_s (Jefferys)
#tauS: prop to 1

#some global parameters
n = dim(math)[1]
y = math$mathscore
p=length(table(math$school))

#MCMC settings
NMC = 21000
burn = 1000

#store sampling results
mu_save = rep(0,NMC)
sigmaS_save = rep(0,NMC)
tauS_save = rep(0,NMC)
theta_save = matrix(0,NMC,p)

#initial parameters for sampling
sigmaS = 1
tauS = 1
mu = mean(math$mathscore)
#mu = 0

#index matrix
M = genInd(math$school,c(1:100))

for (i in 1:(NMC+burn)) {
  if(i %% 250 == 0) cat(paste0("Iteration ", i, "\n"))
  
  m = mu*rep(1,p)
  
  #update theta, given y, mu,sigmaS, tauS
  t.cov = solve((1/sigmaS)*crossprod(M,M)+(tauS*sigmaS)^(-1)*diag(p))
  t.avg = crossprod(t.cov, (crossprod(M, (1/sigmaS)*y)+((tauS*sigmaS)^(-1)*diag(p)) %*% m))
  #t.avg = t(t.cov) %*% (t(M) %*% ((1/sigmaS)*diag(n)) %*% y + ((tauS*sigmaS)^(-1)*diag(p)) %*% m)
  theta = rmvnorm(1, t.avg, t.cov)
  
  all.avg = theta[math$school]
  
  #update mu, given y ,tauS,sigmaS, theta
  mu.var = tauS*sigmaS/p
  mu.avg = sum(theta)/p
  mu.sd = sqrt(mu.var)
  mu = rnorm(1,mean = mu.avg, sd = mu.sd)
  
  #update sigmS, given y, mu,tauS, theta
  sigmaS.a = 0.5*(n+p)
  sigmaS.b = 0.5*(crossprod((y-all.avg),(y-all.avg))+(1/tauS)*sum((theta-m)^2))
  sigmaS = rinvgamma(1,sigmaS.a,sigmaS.b)
  
  #update tauS, given y, sigmaS, mu, theta
  tauS.a = 0.5*p-1
  tauS.b = 0.5*(1/sigmaS)*sum((theta-m)^2)
  tauS = rinvgamma(1,tauS.a, tauS.b)
  
  
  if (i > burn) {
    mu_save[i-burn] = mu
    sigmaS_save[i-burn] = sigmaS
    tauS_save[i-burn] = tauS
    theta_save[i-burn,] = theta
  }
  
}

#analyze shrinkage
t.pred =colMeans(theta_save)
sk=abs((t.pred-y.bar$mathscore)/y.bar$mathscore)
plot(sample.num,sk)



