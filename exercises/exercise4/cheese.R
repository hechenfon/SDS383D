#This is for exercise 4, section 2
library(lme4)
library(MCMCpack)
library(mvtnorm)
library(lattice)
library(ggplot2)

rm(list=ls())


getSampleThetaI <- function(i, group, x, y, vS, sigma, mu) {
   xi = x[group==i,]
   yi = y[group==i]
   rauMat = (1/vS)*diag(dim(xi)[1])
   prec = solve(sigma)
   covMat = solve(crossprod(xi, rauMat %*% xi)+prec)
   mm = t((crossprod(yi, rauMat %*% xi)+crossprod(mu, prec)) %*% covMat)
   ss = t(rmvnorm(1,mm, covMat))
  return(ss)
}

############read in data###############
cheese<-read.csv('cheese.csv')
cheese$logP=log(cheese$price)
cheese$logQ<-log(cheese$vol)


################################################################
########fit linear model with 'lm'##############################
################################################################
#figure 1: with/without display, 
color = c("red","black")
plot(cheese$logP,cheese$logQ,col=color[factor(cheese$disp)],xlab="log of price",ylab="log of sale volume")
par(xpd=TRUE)
legend(1.5,20,legend=c("no-display","with-display"),col=color,pch=21)
par(xpd=FALSE)

#visualize data per store
xyplot(logQ~logP|store,data=cheese, col=c('blue','red'),group=disp)

#complete-pooling
comPlm<-lm(logQ~logP,data=cheese)

#no-pooling fitting
noPlm<-lm(logQ~logP+factor(store)+disp+logP:disp-1,data=cheese)  #linear regression, each store different slope, considering display
plot(cheese$logQ,predict(noPlm))
abline(0,1,col="red")

################################################################
#############fit linear hierachical model with'lmer'############
################################################################
#lmer for hierachical ReML estimation
hlm0<-lmer(logQ~logP+disp+(1|store),data=cheese)

hlm1<-lmer(logQ~logP+disp+(1+logP|store),data=cheese)

hlm2<-lmer(logQ~logP+disp+(1+logP+disp|store),data=cheese)

hlm3<-lmer(logQ~(1+logP+disp+disp:logP|store),data=cheese)

plot(cheese$logQ,predict(hlm3))
abline(0,1,col="red")

#########################################################
#########fit hlm3 with full bayesian way#################
#########################################################
#prior setting:
##vS(random effect): prop to 1
##(mu, sigma): ~NIW(m0, w0, c0, d0)

group = as.numeric(factor(cheese$store))
X = as.matrix(data.frame(cheese$logP, cheese$disp, cheese$disp*cheese$logP, 1))
y = cheese$logQ
idx = as.numeric(names(table(as.numeric(group))))
#idx = names(table(factor(cheese$store)))

p = dim(X)[2]
n = dim(X)[1]
g = length(table(group))

#set prior parameters:
m0 = rep(1,p)*0.001
w0 = 0.001
c0 = diag(p)*0.001
d0 = 1

######################################################
#######################gibbs sampling#################
######################################################

#MCMC settings
NMC = 5000
burn = 1000

#store sampling results
mu_save = matrix(0,NMC,p)
sigma_save = matrix(0, NMC, p^2)
vS_save = rep(0, NMC)
theta_save = matrix(0,NMC,g*p)
rss_save = rep(0,NMC+burn)

#initial parameters for sampling
#mu = colMeans(X)
mu = rep(0,p)
sigma = diag(p)
vS = 1

for (i in 1:(NMC+burn)) {
  if(i %% 250 == 0) cat(paste0("Iteration ", i, "\n"))
  
  #######update theta, given y, mu, sigma, vS######
  #each theta[i] is mvn
  theta = rep(0,g*p)
  theta.idx = sapply(idx,getSampleThetaI,x=X, group=group, y = y, vS=vS, sigma=sigma, mu=mu)
  theta = matrix(theta.idx, 1, g*p)
  
  ######update (mu,sigma) given theta(tM)######
  #Normal Inverse Wishart
  tM = matrix(theta,nrow = g, ncol = p, byrow=TRUE)
  #posterior update P((mu,sigma)|theta)
  #w*
  w2 = w0+g
  #d*
  d2 = d0+g
  #c*
  x.bar = colMeans(tM)
  s = crossprod((tM-x.bar),(tM-x.bar))
  c2 = c0+s+((w0*g)/(w0+g))*((x.bar-m0) %*% t(x.bar-m0))
  #m*
  m2 = (w0*m0+g*x.bar)/(w0+g)
  #sample from posterior distribution NIW
  sigma = riwish(d2,c2)
  mu = t(rmvnorm(1,m2,(1/w2)*sigma))
  
  ######update vS: Inverse Gamma#######
  thetaM = tM[group,]
  #rss = sum((y-diag(X %*% t(thetaM)))^2) #slow way
  rss = sum((y-rowSums(X * thetaM))^2)
  ivs = rgamma(1, (0.5*n-1), 0.5*rss)
  vS = 1/ivs
  
  if (i > burn) {
    theta_save[i-burn,] = theta
    mu_save[i-burn,] = t(mu)
    sigma_save[i-burn,] = matrix(sigma,1,p^2)
    vS_save[i-burn] = vS
    
  }
  rss_save[i] = rss
}

######posterior distribution analysis based on sampling######
theta.m = colMeans(theta_save)
theta.m = matrix(theta.m, nrow = g, ncol = p, byrow = TRUE)
y.pred = rowSums(X * theta.m[group,])

y.bayes.resid = abs(y-y.pred)
y.hlm.resid = abs(resid(hlm3))
y.noPlm.resid = abs(resid(noPlm))

######compare the residue with bayesian/lmer()######
plot(y.bayes.resid,y.hlm.resid)
#plot(y.bayes.resid,y.noPlm.resid)
abline(0,1,col="red")

theta.p = colMeans(theta)
#tM.p = as.data.frame(matrix(theta.p,nrow = g, ncol = p, byrow=TRUE))
tM.p = as.data.frame(coef(hlm3)$store)
colnames(tM.p) = c('logP','disp','disp.logP','intercept')
tM.p$store = idx
#tM.p$disp = factor(tM.p$disp)

scplot<-ggplot(cheese, aes(logP, logQ))+
  geom_point(pch=1,aes(colour=factor(disp)))+
  facet_wrap(~factor(as.numeric(factor(store))),ncol=9)+
  scale_colour_manual("Display", values = c("blue", "red"))

abplot<-scplot+geom_abline(data=tM.p, aes(slope=logP,intercept=intercept),col='blue')+
    geom_abline(data=tM.p, aes(slope=logP+disp.logP,intercept=intercept),col='red')

pdf('post.pred1.pdf',width=20,height=12)
abplot
dev.off()
