myPalette <- c("#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF")
palette(myPalette); par(pch=19)
plotIm <- function(x){image(t(matrix(t(x),ncol=16,byrow=TRUE)[16:1,]),col=gray(255:0/255),axes=F); box()}

###########################################################################
# Loading of data 
# load('data/Wine.Rdata')           # n = 150,  p = 13
# load('data/Chironomus.Rdata')     # n = 149, p = 17
# load('data/USPS358.Rdata')        # n = 1756, p = 256
# load('data/Mars.Rdata')           # n = 38400, p = 255
# load('data/Galaxy-small.Rdata')   # n = 10000, p = 1539
# load('data/NIR_data.Rdata')       # n = 202, p = 2800
# also Velib data

###########################################################################
# Hyper-sphere
p = 1:100
plot(p,pi^(p/2)/gamma(p/2+1),type='b',col=1)
title(main="Volume of unit hyper-sphere")

###########################################################################
#Shell
plot(p,1-0.9^p,type='b',col=1)
title(main="P(X in C_0.9)")

############################################################################
# Bayes classifier vs EDDA
library(mvtnorm); library(MASS); library(Rmixmod)

BayesClass <- function(X,tau,m,S){
  G = length(m); d = ncol(X)
  P = matrix(NA,nrow(X),G)
  for (g in 1:G) P[,g] = tau[g] * dmvnorm(X,m[[g]],S[[g]])
  P = P / rowSums(P) %*% matrix(1,1,G)
  list(P = P, cls = max.col(P))
}

n = 120; nbrep = 5
dims = seq(10,210,20)
err = err2 = err3 = matrix(NA,length(dims),nbrep)
for (i in 1:length(dims)){
  cat('.')
  for (j in 1:nbrep){
    # Simulation
    d = dims[i]
    m1 = c(0,0,rep(0,d-2)); m2 = c(0,-2,rep(0,d-2)); m3 = c(2,0,rep(0,d-2));
    S1 = diag(d); S2 = 0.8 * diag(d); S3 = 0.9 * diag(d)
    X = as.data.frame(rbind(mvrnorm(n/3,m1,S1),mvrnorm(n/3,m2,S2),mvrnorm(n/3,m3,S3)))
    X2 = as.data.frame(rbind(mvrnorm(10*n/3,m1,S1),mvrnorm(10*n/3,m2,S2),mvrnorm(10*n/3,m3,S3)))
    cls = rep(1:3,rep(n/3,3))  
    cls2 = rep(1:3,rep(10*n/3,3))
    # Classification with the Bayes' classifier
    pred = BayesClass(X2,rep(1/3,3),list(m1,m2,m3),list(S1,S2,S3))
    # Classification with EDDA
    mod = mixmodLearn(X,cls,models=mixmodGaussianModel(listModels = 'Gaussian_pk_Lk_I'))
    res = mixmodPredict(X2,mod["bestResult"])@partition
    # Computing error rate
    err[i,j] = sum(pred$cls != cls2) / length(cls2)
    err2[i,j] = sum(res != cls2) / length(cls2)
  }
}
cat('\n')
boxplot(t(err),ylim=c(0.1,0.33),names=dims,xlab='Dimension',ylab='Classification error rate',col=3)
boxplot(t(err2),names=dims,xlab='Dimension',ylab='Classification error rate',col=4,add=TRUE)
legend("bottomleft",legend = c('Bayes classifier','EDDA'),col=c(3,4),lty=1,pch=19)

###########################################################################
# PCA basics
load('data/Wine.Rdata')
pc = princomp(X)
plot(pc)
biplot(pc)
pairs(predict(pc)[,1:5],col=cls,pch=19)

###########################################################################
# PCA on NIR data
load('data/NIR_data.Rdata')
matplot(t(Y),col=1,type='l',lty=(3:1)[cls],xlab='Wavelengths',ylab='Intensity')
matplot(t(Y),col=cls,type='l',lty=(3:1)[cls],xlab='Wavelengths',ylab='Intensity')

#pc = princomp(Y)

library(MASS)
U = svd(Y)$v
X = as.matrix(Y) %*% U
pairs(X[,1:4],col=cls,pch=(15:17)[cls],cex=1.25,labels=c('PC 1','PC 2','PC 3','PC4'))

###########################################################################
# PCA vs FDA on USPS data
load('data/USPS358.Rdata')

pc = princomp(X)
plot(predict(pc),col=cls)

out = lda(X,cls)
Xproj = as.matrix(X) %*% out$scaling
plot(Xproj,col=cls)

###########################################################################
# GMM with Mclust on Wine and Chironomus
# Try also on NIR data or USPS
library(mclust)
load('data/Wine.Rdata')
out = Mclust(X,G=3)
table(out$classification,cls)
adjustedRandIndex(out$classification,cls)
plot(predict(princomp(X)),col=out$class)

###########################################################################
# GMM with Mclust on principal components
Xp = predict(princomp(X))[,1:2]
out = Mclust(Xp,G=3)
table(out$classification,cls)
adjustedRandIndex(out$classification,cls)
plot(predict(princomp(X)),col=out$class)

###########################################################################
# HDDC (look at parameters and choice of the dimensionality)
library(HDclassif)
load('data/Wine.Rdata')
out = hddc(X,3)
adjustedRandIndex(out$class,cls)
plot(predict(princomp(X)),col=out$class)

###########################################################################
# FisherEM (look at parameters and visualization)
library(FisherEM)
out = fem(X,3)
adjustedRandIndex(out$cls,cls)
plot(predict(princomp(X)),col=out$cls)
plot(out)

###########################################################################
# SparseFEM
out = sfem(X,3,l1 = 0.3,model='all')
adjustedRandIndex(out$cls,cls)
plot(out)

class(out$U) = "loadings"
out$U
