###########################################################
########## Created by Chenfeng He on 1 May 2017 ###########
###########################################################

#Data resource: Guo G, Huss M, Tong G Q, et al. Resolution of cell fate decisions revealed by single-cell gene expression analysis from zygote to blastocyst[J]. Developmental cell, 2010, 18(4): 675-685.
#Data description: 48 genes expression value from (1->2->4->8->16->32->64 cells)
#Plan: map high-dimensional cells into lower dimensional space with:
# 1. PCA
# 2. FA(factor analysis) with public R package
# 3. FA with bayesian method(Rowe D B, Press S J. Bayesian Factor Analysis By Gibbs Sampling and Iterated Conditional Modes[J]. 1994.)
library(ggfortify)
library(psych)
library(covreg)
library(ggplot2)
library(MCMCpack)

#house clean
rm(list=ls())

#read in data
gexp <- read.csv('./data/mmc4.csv', header=T, row.names = 1, stringsAsFactors = F)

#normalize by subtracting control genes(average of Actb and Gapdh )
gexp <- gexp-(gexp$Actb+gexp$Gapdh)/2
#remove control genes
gexp <- gexp[,colnames(gexp) != 'Gapdh' & colnames(gexp) != 'Actb']
#normalize
gexp <- scale(gexp,center=T,scale=T)

#64 cell data
g64.exp <- gexp[grepl('64C',row.names(gexp)),]

#heatmap on 64 cell level
hm<-heatmap(as.matrix(g64.exp))

#add biological label(64 cells based on hierarchical clustering, 1...32 based on experiment)
cell.grp<-as.data.frame(cutree(hclust(dist(g64.exp)),k=3))

gexp<-merge(gexp,cell.grp,by=0,all=TRUE)

gexp[grepl('^1C',gexp[,1]),48]<-'1.c'
gexp[grepl('^2C',gexp[,1]),48]<-'2.c'
gexp[grepl('^4C',gexp[,1]),48]<-'4.c'
gexp[grepl('^8C',gexp[,1]),48]<-'8.c'
gexp[grepl('^16C',gexp[,1]),48]<-'16.c'
gexp[grepl('^32C',gexp[,1]),48]<-'32.c'
gexp[gexp[,48]==1,48]<-'TE'
gexp[gexp[,48]==2,48]<-'PE'
gexp[gexp[,48]==3,48]<-'EPI'

names(gexp)[48] <- 'bio.lab'
rownames(gexp) = gexp[,1]
gexp<-gexp[,-c(1)]

#PCA analysis for 64 cells level
g64.exp <- gexp[grepl('64C',row.names(gexp)),]
g64.pca<-prcomp(g64.exp[,1:46]) 
summary(g64.pca)  #first 2 pc variance: 0.5347+0.1327
autoplot(g64.pca,data=g64.exp,colour='bio.lab')
#PCA map all cells
gexp.rv = as.matrix(gexp[,1:46]) %*% g64.pca$rotation
rv.lab = data.frame(gexp.rv[,1:2],gexp$bio.lab)
colnames(rv.lab) = c('PC1','PC2','bio.lab')
#map all cells on the PC1-PC2 subspace
#qplot(PC1,PC2,data=rv.lab,colour=bio.lab)
ggplot(rv.lab,aes(x=PC1,y=PC2,color=bio.lab))+geom_point()+theme_bw()


###public available FA package for factor analysis###
fld<-fa(g64.exp[,1:46],nfactors = 2,rotate='varimax',fm="minres")
fld.res<-data.frame(fld$scores,g64.exp$bio.lab)
#ld <- fld$loadings
#fld.res<-data.frame(t(solve(t(ld) %*% ld) %*% t(ld) %*% t(as.matrix(gexp[,1:46]))),gexp$bio.lab)
colnames(fld.res)<-c('FA1','FA2','bio.lab')
ggplot(fld.res,aes(x=FA1,y=FA2,color=bio.lab))+geom_point()+theme_bw()




#####GIBBS start here######
X = as.matrix(g64.exp[,1:46])
n = dim(X)[1]   #number of observation
p = dim(X)[2]   #number of original dimension
m = 2           #number of target dimension

###set up parameters###
NMC = 50000
burn = 30000
thin = 2

#save sampled values
lb.save = array(1,dim=c(p,m,NMC+burn))
phi.save = array(1,dim=c(p,p,NMC+burn))
f.save = array(1,dim=c(n,m,NMC+burn))

#hyperparameters for priors
v = 30
B = 0.2*diag(p)
lb0 = fld$loadings
H = 10*diag(m)

#initial values
lb.save[,,1] = fld$loadings
phi.save[,,1] = diag(p)
f.save[,,1] = fld$scores

#sampling
for (iter in 2:(NMC+burn)) {
  if(iter %% 200 == 0) cat(paste0("Iteration ", iter, "\n"))
  
  #sample on cond. post. for lambda - matrix variate normal 
  tc = solve(H+crossprod(f.save[,,iter-1]))
  m.lb = (t(X) %*% f.save[,,iter-1] + lb0 %*% H) %*% tc
  lb.save[,,iter] = rmn(M = m.lb, Srow = phi.save[,,iter-1], Scol = tc)
  
  #sample on cond. post. phi~ IW
  tc1 = X-f.save[,,iter-1] %*% t(lb.save[,,iter-1])
  tc2 = lb.save[,,iter-1]-lb0
  u = t(tc1) %*% (tc1)+(tc2) %*% H %*% t(tc2)+B
  #u = crossprod(tc1)
  phi.save[,,iter] = riwish(n+m+v-p-1, u)
  
  #sample on cond. post. F ~ matrix variate normal
  tc = solve(diag(m) + t(lb.save[,,iter-1]) %*% solve(phi.save[,,iter-1]) %*% lb.save[,,iter-1])
  m.f = X %*% solve(phi.save[,,iter-1]) %*% lb.save[,,iter-1] %*% tc
  f.save[,,iter] = rmn(M = m.f, Scol = tc, Srow = diag(n))
  
}

#thin the chain
sel.id = seq(burn,NMC+burn,thin)

#analyze gibbs results
#brief check mix
plot(lb.save[1,1,sel.id])

#calculate posterior means
lb.pm = apply(lb.save[,,sel.id],c(1,2),mean)
phi.pm = apply(phi.save[,,sel.id],c(1,2),mean)
f.pm = apply(f.save[,,sel.id],c(1,2),mean)

dt<-data.frame(f.pm,g64.exp$bio.lab)
names(dt)<-c('FA1','FA2','bio.lab')
ggplot(dt,aes(x=FA1,y=FA2,color=bio.lab))+geom_point()+theme_bw()


#varimax rotation