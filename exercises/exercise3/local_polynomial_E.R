#for Excercise 3, Local polynomial regression section, problem E

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

LOOCV.bdw <-function(x,y,h) {
  #leave one out CV, in order to choose best bandwidth
  #in order to minimize rss, gaussian kernel, local linear smoothing
  
  #calculating smoothing matrix first
  len = length(y)
  sm.mat <-matrix(0,len,len)
  for (i in c(1:len)) {
    x.test = x[i]
    kns = sapply((x.test-x)/h,dnorm)
    s1 = kns %*% (x-x.test)
    s2 = kns %*% ((x-x.test)^2)
    w = diag(kns %*% t(s2-s1*(x-x.test)))
    
    sm.mat[i,] = w/sum(w)
  }

  y.pred = sm.mat %*% y
  err = sum(((y-y.pred)*(1/(1-diag(sm.mat))))^2)
  
  return(err)
}


#use data from 'utilities.csv' for test
data<-read.csv("utilities.csv",header=T)
y = data$gasbill/data$billingdays
x = data$temp

######PART I: LOOCV for optimized bandwidth######
# #LOOCV to determine the best bandwidth
# h = c(10:200)/10
# err = sapply(h,LOOCV.bdw, x=x, y=y)
# plot(h,err,xlab='bandwidth',ylab='LOOCV error')

######PART II: after chosen bandwidth, Local linear smooth######
h.opt = 6.9
x.test <- c(10:80)

y.hat = sapply(x,localLinearSmoother,x=x,y=y,h=h.opt)

#plot figures
plot(x,y,xlab='average temperature',ylab='average daily gas bill')
lines(x.test,y.hat)

