# copyright belongs to Sungkyu Jung
# See paper "Kurtosis test of modality for rotationally symmetric distributions on hyperspheres"
# converted code from MATLAB
# output is the type of fitted circle
# likelihood ratio test from shapes package

LRTpval2 <- function(resGreat,resSmall,n) {
  chi2 <- max(n*log(sum(resGreat^2)/sum(resSmall^2)))
  pval <- 1-pchisq(q = chi2, df = 1, lower.tail = T) # likelihood test p-value Also you can use chi2cdf(chi2,1) from library(PEIP) like matlab
}

kurtosisTestFunction <- function(sphericalData, alpha=0.1) {
  ndata<-dim(sphericalData)[2]
  
  subsphereSmall<-getSubSphere(sphericalData,geodesic = "small")
  subsphereGreat<-getSubSphere(sphericalData,geodesic = "great")
  
  currentSphere<-sphericalData
  
  rSmall<-subsphereSmall$r                    
  centerSmall<-subsphereSmall$center      # NB! center is the centerSmall is pnsSmall$PNS$orthaxis[[1]]
  # and centers in matlab
  resSmall <- acos(t(centerSmall)%*%currentSphere)-rSmall  # NB!!! resSmall==(pnsSmall$resmat)[2,] i.e., residuals are second coordinates of PNS
  
  rGreat<-subsphereGreat$r                    
  centerGreat<-subsphereGreat$center          
  resGreat <- acos(t(centerGreat)%*%currentSphere)-rGreat  # NB!!! resGreat==(pnsGreat$resmat)[2,] i.e., residuals are second coordinates of PNS
  
  # LRTpval is the likelihood ratio test from 'shapes' package
  # Chi-squared statistic for a likelihood test
  pval1 <- LRTpval(resGreat,resSmall,n = ndata)
  pval1
  
  if(pval1>alpha){
    print('great by likelihood ratio test')
    return('great')
    break
  }
  
  # # equivalently we can find pval by pns function
  # pnsTest2<-pns(sphericalData)
  # pnsTest2$PNS$pvalues
  # sum(pnsTest2$resmat[2,]==resSmall)
  
  # kurtosis test routine
  X <- LogNPd(rotMat(centerSmall) %*% currentSphere)
  
  # Note that the tangential point is the center of the small circle
  d<-dim(X)[1]
  n<-dim(X)[2]
  normX2 <- colSums(X^2)
  kurtosis <- sum( normX2^2 ) / n / ( sum( normX2 ) / (d * (n-1)) )^2
  M_kurt <- d * (d+2)^2 / (d+4)
  V_kurt <- (1/n) * (128*d*(d+2)^4) / ((d+4)^3*(d+6)*(d+8))
  pval2 <- pnorm((kurtosis - M_kurt) / sqrt(V_kurt))
  
  if(pval2>alpha){
    return('great')
  }else{
    # drawCircleS2(normalVec = centerSmall,radius = rSmall)
    return('small')
  }
}