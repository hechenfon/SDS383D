#this is for Excercise 3, local polynomial section, problem F
rm(list=ls())


localLinearSmoother <- function(x.test,h) {
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
x.test = x
y.hat1 = sapply(x.test,localLinearSmoother,h=h.opt)
residue = y.hat1-y
rss1 = sum((y.hat1-y)^2)

tmp<-data.frame(x,residue)
tmp$bins<-cut(tmp$x,breaks=c(0:10)*10,labels=c(1:10))
cuts<-sapply(c(1:10),function (b) {
  a = tmp[tmp$bins==b,2]
  return(sd(a))
}
)
cuts[is.na(cuts)]=min(cuts[!is.na(cuts)])
h=cuts/min(cuts)
tmp$bdw = sapply(tmp$bins, function(b){tmp$bdw[tmp$bins==b] = h[b] })
y.hat2=mapply(localLinearSmoother,x.test=tmp$x,h=tmp$bdw)
rss2 = sum((y.hat2-y)^2)



