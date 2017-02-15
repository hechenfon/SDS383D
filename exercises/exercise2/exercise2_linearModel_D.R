#this is for exercise 2, conjugate gaussian linear model, basic section, problem D

gdp<-read.csv("gdpgrowth.csv",header=T,stringsAsFactors = FALSE)

X = cbind(rep(1,dim(gdp)[1]),gdp$DEF60) #add one column for interception
y = gdp$GR6096

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
m2 = solve(k2) %*% (t(X) %*% y+k %*% m)
#eta*
eta2 = eta+t(y) %*% y+t(m) %*% k %*% m-t(m2) %*% k2 %*% m2

#posterior distribution of beta: p(beta|y)
#according to (C), p(beta|y) is a t-distribution: t(m2,d2,(eta2/d2)K2), posterior estimate of beta is m2
beta_post= m2

pred_y = X %*% beta_post

