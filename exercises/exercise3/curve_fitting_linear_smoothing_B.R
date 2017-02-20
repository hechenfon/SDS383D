#this is for Exercise 3, 'Curve fitting by linear smoothing', problem B

norm.kernel.sm<-function(dist, bdw) {
  #normal kernel function
  #dist: distance
  #bdw: bandwidth
  rt = dist/bdw
  #for normal kernel
  pos = (1/sqrt(2*pi))*exp(-0.5*rt^2)
  return(pos)
}

wgt.fun<-function(x,x.new,myfun,bdw) {
  #calculate weight parameters
  #input
    #x: predictor matrix
    #x.new: new sample
    #myfun: name of kernel smooth function
    #bdw: smoother bandwidth
  #output
    #a vector contains all the weights
  kernel.sm = match.fun(myfun)
  n = dim(x)[1]
  p = dim(x)[2]
  dist = sqrt(rowSums((x - t(x.new %*% matrix(rep(1,n),ncol=n)))^2))
  wt = sapply(dist,kernel.sm,bdw=bdw)
  #normalize weight parameters
  wt = wt/sum(wt)
  return(wt)
}

#set up test samples, test on one dimension , f(x)=sin(x), add gaussian noise
x = matrix(runif(100,min=-10,max=10),nrow = 100)
y = sin(x)+rnorm(100,mean=0,sd=0.5)
  #bandwidth
h = 1

xx = seq(-10,10,0.1)
test.x = matrix(xx,nrow=length(xx))

beta=sapply(test.x,wgt.fun,x = x, myfun=norm.kernel.sm,bdw=h)
y2 = t(y) %*% beta

