myPalette <- c("black","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF")
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
# Hyper-sphere and shell
p = 1:100
plot(p,pi^(p/2)/gamma(p/2+1),type='b',main="Volume of the hyper-sphere")

plot(p,1-0.9^p,type='b',main='P(X in S0.9)')
###########################################################################
# PCA basics
load('data/Wine.Rdata') 
pc  = princomp(X)

S = cov(X)
vect = eigen(S)$vect

Xp = predict(pc)[,1:2]
Xp = as.matrix(X) %*% vect[,1:2]
plot(Xp,col=cls)
###########################################################################
# PCA on NIR data
load('data/NIR_data.Rdata')
matplot(t(Y),col=cls,type='l',lty=1)

#pc = princomp(Y) # does not work since n < p
SS = cov(t(Y)) # PCA through eigen(X %*% t(X))
V = eigen(SS)$vect
U = t(Y) %*% V
Yp = as.matrix(Y) %*% U[,1:10]
###########################################################################
# PCA vs FDA on USPS data
load('data/Wine.Rdata') 
library(MASS)
out = lda(X,cls)
Xp = as.matrix(X) %*% out$scaling
plot(Xp,col=cls)

###########################################################################
# GMM with Mclust on Wine and Chironomus
# Try also on NIR data or USPS
library(mclust)
load('data/Wine.Rdata')
out = Mclust(X,3,modelNames = 'VVV')
adjustedRandIndex(out$class,cls)

###########################################################################
# GMM with Mclust on principal components (HOMEWORK)


###########################################################################
# HDDC (look at parameters and choice of the dimensionality)
#install.packages('HDclassif')
library(HDclassif)
load('data/Wine.Rdata')
?hddc
out = hddc(X,3)

Q = out$Q
par(mfrow=c(1,3))
for (k in 1:3){
  Qk = Q[[k]]
  Xp = as.matrix(X) %*% Qk[,1:2]
  plot(Xp,col=cls)
}

library(mclust)
adjustedRandIndex(out$class,cls)

out = hddc(X,3,model='all')
adjustedRandIndex(out$class,cls)

out = hddc(X,K=2:7)
adjustedRandIndex(out$class,cls)
###########################################################################
# FisherEM (look at parameters and visualization)
#install.packages('FisherEM')
library(FisherEM)
?fem
out = fem(X,3)
adjustedRandIndex(out$cls,cls)

out = fem(X,3,model='all')
adjustedRandIndex(out$cls,cls)

plot(out)
###########################################################################
# SparseFEM
out = sfem(X,3,l1 = 0.4,model='all')
adjustedRandIndex(out$cls,cls)
out$U
plot(out)
