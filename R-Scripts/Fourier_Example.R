#################################################################
#################################################################
# R example

library(fda)
library(refund)
library(ggplot2)
library(dplyr)
library(reshape2)

# basis functions
fourier.basis<-create.fourier.basis(rangeval=c(0,10), nbasis=5)
plot(fourier.basis, lty=1, lwd=2)
par(mfrow=c(1,1))

# Example based on Brownian motion
set.seed(1234)
Nsim<-1
Times<-10000
W.mat=matrix(0, ncol=Nsim, nrow=Times)
for(n in 1:Nsim){W.mat[, n]=cumsum(rnorm(Times))/100}

plot.new()
for (i in 1:Nsim) {
  plot(1:Times,W.mat[,i],type = 'l',
       col=sample(1:Times,size = 1),ylim = c(min(W.mat),max(W.mat)),xlab = '',ylab = '') 
  par(new=TRUE)
}

plot.new()

# smooth by fourier basis
f.basis=create.fourier.basis(rangeval=c(0,Times), nbasis=40)
W.fd=smooth.basis(y=W.mat, fdParobj=f.basis)

plot(W.fd$fd, ylab="", xlab="",col='blue',lty=1)



# Mean in case of Nsim > 1
if(Nsim>1){
  W.mean=mean.fd(W.fd$fd)
  lines(W.mean, lty=2, lwd=3) 
}
