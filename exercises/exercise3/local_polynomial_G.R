#this is for Excercise 3, local polynomial section, problem G, 
#95% confidence interval

rm(list=ls())
localLinearSmoother <- function(x.test,x,y,h) {
  #This function is to do local linear estimate, 
  #given observed values and bandwidth, 
  #use gaussian kernel and return predicted values
  
  #given x, get smooth kernel density
  kns = sapply((x.test-x)/h,dnorm)
  s1 = kns %*% (x-x.test)
  s2 = kns %*% ((x-x.test)^2)
  w = diag(kns %*% t(s2-s1*(x-x.test)))
  
  y.pred = (w %*% y)/sum(w)
  
  return(y.pred)
}

#use data from 'utilities.csv' for test
data<-read.csv("utilities.csv",header=T)
y = data$gasbill/data$billingdays
x = data$temp

h.opt = 6.9
x.test <- c(10:80)
n = length(y)

y.hat = sapply(x.test,localLinearSmoother,x=x,y=y,h=h.opt)
y.low95CI = y.hat-1.96*(h.opt/sqrt(n))
y.up95CI = y.hat+1.96*(h.opt/sqrt(n))
#plot figures
plot(x,y,xlab='average temperature',ylab='average daily gas bill')
lines(x.test,y.hat)
lines(x.test,y.low95CI,lty=2,col='red')
lines(x.test,y.up95CI,lty=2,col='red')



