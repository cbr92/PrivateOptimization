#############################################
########## CONTENTS #########################
#############################################

## 1. implementations of huber loss and its first derivative
## 2. implementations of self-concordant loss functions from 
##    Ostrovskii & Bach (2021)** and their derivatives
## 3. generate Figure [4]

#**Ostrovskii, D. M., & Bach, F. (2021). Finite-sample analysis of M-estimators using self-concordance. 
#Electronic Journal of Statistics, 15(1), 326-391.

#############################################
##################   1   ####################
#############################################
huberloss<-function(c,a){
  if(abs(a)<=c){
    return((a^2)/2)
  }
  else{
    return(c*abs(a)-(c^2)/2)
  }
}

psiHuber <- function (u, k) 
{
  u*pmin(1, k/abs(u))
} 


#############################################
##################   2   ####################
#############################################

# a (2c/3)-self concordant loss function from
# Ostrovskii & Bach (2021), equation (21)
SCloss1<-function(v,c=1){
  v<-v/c
  r<-sqrt(1+4*v^2)
  calc<-(r-1+log((r-1)/(2*v^2)))/2
  calc[v==0]<-0
  
  return((c^2)*calc)
}

SClossderiv1<-function(v,c=1){
  v<-v/c
  r<-sqrt(1+4*v^2)
  calc<-2*v/(r-1)-1/v
  calc[v==0]<-0
  calc[abs(v)<0.00002]<-0
  return(c*calc)
}

# a (3c/2)-self concordant loss function from
# Ostrovskii & Bach (2021), equation (19)
SCloss2<-function(v,c=1){
  v<-v/c
  calc<-log(cosh(v))
  return((c^2)*calc)
}

SCloss2deriv<-function(v,c=1){
  v<-v/c
  calc<-tanh(v)
  return(c*calc)
}

#############################################
##################   3   ####################
#############################################
source("Users/caseybradshaw/Documents/privacy/SmoothedHuber.R")

x<-seq(from=-2,to=2,by=0.01)

huber<-sapply(x,huberloss,c=1.345)
smoothed.huber<-sm.Huber(x,k=1.345,h=0.5)
SChuber1<-SCloss1(v=x,c=1.345)
SChuber2<-SCloss2(v=x,c=1.345)

huberderiv<-sapply(x,psiHuber,k=1.345)
smoothed.huberderiv<-sm.psi(x,k=1.345,h=0.5)
SChuber1deriv<-SClossderiv1(v=x,c=1.345)
SChuber2deriv<-SCloss2deriv(v=x,c=1.345)

par(mfrow=c(1,2))
plot(x,huber,type="l",xlab="x",ylab="Loss",ylim=c(0,2),main="Example Loss Functions",lwd=1.5)
lines(x,smoothed.huber,col=2,lwd=1.5)
lines(x,SChuber1,col=3,lwd=1.5)
lines(x,SChuber2,col=4,lwd=1.5)
#legend("top",legend=c("Smoothed Huber loss","Huber loss",expression(paste(c^2,"log(cosh(x/c))")),"(2/c,3)-self concordant"), 
#       lty=1,col=c(2,1,4,3),bty="n",cex=0.8)
legend("top",legend=c("Smoothed Huber loss","Huber loss","(3/c,2)-self concordant","(2/c,3)-self concordant"), 
       lty=1,lwd=1.5,col=c(2,1,4,3),bty="n",cex=0.8)


plot(x,huberderiv,xlab="x",ylim=c(-2,2),main="Loss Function Derivatives",type="l",ylab="",lwd=1.5)
lines(x,smoothed.huberderiv,col=2,lwd=1.5)
lines(x,SChuber1deriv,col=3,lwd=1.5)
lines(x,SChuber2deriv,col=4,lwd=1.5)
legend("topleft",legend=c("Smoothed Huber loss","Huber loss","(3/c,2)-self concordant","(2/c,3)-self concordant"), 
       lty=1,lwd=1.5,col=c(2,1,4,3),bty="n",cex=0.8)
