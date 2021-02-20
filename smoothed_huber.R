#### Functions for a "smoothed" version of Huber loss
#### specifically motivated by wanting the second
#### derivative of the Huber loss decay smoothly 
#### from 1 to zero over an interval of length "h".
#### This is achieved by using a cubic polynomial
#### to interpolate between the points (k,1)
#### and (k+h,0) (and analogously for neative x).
#### This cubic polynomial is unique given the
#### constraint that it be continuously
#### differentiable. From there, the first
#### derivative of the smoothed Huber loss and
#### the smoothed Huber loss itself are also
#### uniquely determined by the requirement that
#### both be continuously differentiable.

#### The functions begin with the second derivative
#### of the smoothed Huber loss, followed by the
#### first derivative of smoothed Huber loss,
#### then the smoothed Huber loss itself. The third
#### derivative is also included for completeness.





#3rd degree polynomial that smooths the corners of the second derivative of Huber loss
P3<-function(x,k,h){
  a<-2/(h^3)
  R<-k-h/2
  
  px<-a*(x-R)*(x-(k+h))^2
  
  return(px)
}

#piecewise function for second derivative of smoothed Huber loss
psiprimeHuber.sm<-function(u,k,h){
  if(abs(u)<k){return(1)}
  else if(abs(u)>(k+h)){return(0)}
  else{
    
    return(P3(x=abs(u),k=k,h=h))
  }
  
}

#for using above function with vector input:
sm.psiprime<-function(x,k,h){
  return(sapply(x,psiprimeHuber.sm,k=k,h=h))
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




#second degree polynomial for interpolated section of third derivative of smoothed Huber loss
P2<-function(x,k,h){
  a<-2/(h^3)
  R<-k-h/2
  
  px<-a*(x-(k+h))^2+2*a*(x-R)*(x-(k+h))
  
  return(px)
}

#piecewise function for third derivative of smoothed Huber loss
psidprimeHuber.sm<-function(u,k,h){
  if(abs(u)<k){return(0)}
  else if(abs(u)>(k+h)){return(0)}
  else{
    return(sign(u)*P2(x=abs(u),k=k,h=h))
  }
}

#for using above function with vector input:
sm.psidprime<-function(x,k,h){
  return(sapply(x,psidprimeHuber.sm,k=k,h=h))
}

