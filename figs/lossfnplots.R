


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


#4th degree polynomial that smooths the corners of the first derivative of Huber loss
P4<-function(x,k,h){
  a<-2/(h^3) #leading coefficient of the cubic polynomial
  R<-k-h/2 #other root of the cubic polynomial (k+h is a root of order 2)
  b<-k+h
  
  d4<-(1/4)*((x^4)-k^4)
  d3<-(-1/3)*((x^3)-k^3)*(R+2*b)
  d2<-(1/2)*((x^2)-k^2)*(((b)^2)+2*R*b)
  d1<-(-1)*(x-k)*R*(b)^2
  
  px<-a*(d4+d3+d2+d1)+k
  
  return(px)
  
}

#piecewise function for derivative of smoothed Huber loss
psiHuber.smoothed<-function(u,k,h){
  if(abs(u)<k){return(u)}
  else if(abs(u)>(k+h)){
    return(sign(u)*P4(x=(k+h),k=k,h=h))
  }
  else{
    return(sign(u)*P4(x=abs(u),k=k,h=h))
  }
}

#for using above function with vector input:
sm.psi<-function(x,k,h){
  return(sapply(x,psiHuber.smoothed,k=k,h=h))
}



# 5th degree polynomial for interpolated section of smoothed Huber loss
P5<-function(x,k,h){
  a<-2/(h^3) #leading coefficient of the cubic polynomial
  R<-k-h/2 #other root of the cubic polynomial (k+h is a root of order 2)
  b<-k+h
  
  d5<-(1/20)*((x^5)-k^5)
  d4<-(-1/12)*((x^4)-k^4)*(R+2*(b))
  d3<-(1/6)*((x^3)-k^3)*(((b)^2)+2*R*(b))
  d2<-(-1/2)*((x^2)-k^2)*R*(b)^2
  d1<-(x-k)*(k-a*((k^4)/4+((k^3)/3)*(-R-2*k-2*h)+((k^2)/2)*((b)*(b+2*R))+k*(-R*(b)^2)))
  
  px<-a*(d5+d4+d3+d2)+d1+(k^2)/2
  
  return(px)
  
}

#piecewise function for smoothed Huber loss
Huber.smoothed<-function(u,k,h){
  if(abs(u)<k){return((u^2)/2)}
  else if(abs(u)>(k+h)){
    return(P5(x=(k+h),k=k,h=h)+P4(x=(k+h),k=k,h=h)*(abs(u)-(k+h)))
  }
  else{
    return(P5(x=abs(u),k=k,h=h))
  }
}

#for using above function with vector input:
sm.Huber<-function(x,k,h){
  return(sapply(x,Huber.smoothed,k=k,h=h))
}





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




x<-seq(from=-2,to=2,by=0.01)

huber<-sapply(x,huberloss,c=1.345)
smoothed.huber<-sm.Huber(x,k=1.345,h=0.5)
SChuber1<-SCloss1(v=x,c=1.345)
SChuber2<-SCloss2(v=x,c=1.345)


plot(x,huber,type="l",xlab="x",ylab="Loss",ylim=c(0,2),main="Example Loss Functions",lwd=1.5)
lines(x,smoothed.huber,col=2,lwd=1.5)
lines(x,SChuber1,col=3,lwd=1.5)
lines(x,SChuber2,col=4,lwd=1.5)
#legend("top",legend=c("Smoothed Huber loss","Huber loss",expression(paste(c^2,"log(cosh(x/c))")),"(2/c,3)-self concordant"), 
#       lty=1,col=c(2,1,4,3),bty="n",cex=0.8)
legend("top",legend=c("Smoothed Huber loss","Huber loss","(3/c,2)-self concordant","(2/c,3)-self concordant"), 
       lty=1,lwd=1.5,col=c(2,1,4,3),bty="n",cex=0.8)


x<-seq(from=-2,to=2,by=0.01)

huberderiv<-sapply(x,psiHuber,k=1.345)
smoothed.huberderiv<-sm.psi(x,k=1.345,h=0.5)
SChuber1deriv<-SClossderiv1(v=x,c=1.345)
SChuber2deriv<-SCloss2deriv(v=x,c=1.345)


plot(x,huberderiv,xlab="x",ylim=c(-2,2),main="Loss Function Derivatives",type="l",ylab="",lwd=1.5)
lines(x,smoothed.huberderiv,col=2,lwd=1.5)
lines(x,SChuber1deriv,col=3,lwd=1.5)
lines(x,SChuber2deriv,col=4,lwd=1.5)
legend("topleft",legend=c("Smoothed Huber loss","Huber loss","(3/c,2)-self concordant","(2/c,3)-self concordant"), 
       lty=1,lwd=1.5,col=c(2,1,4,3),bty="n",cex=0.8)

par(mfrow=c(1,2))
