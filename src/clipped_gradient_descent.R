# contents: 
# noisy gradient descent functions for
# linear regression and logistic regression,
# modified to utilize gradient clipping


#############################################
########### auxiliary functions #############
#############################################
psiHuber <- function (u, k) 
{
  u*pmin(1, k/abs(u))
} 


Fisher.constant <- function(k){ # computes Fisher consistency constant for Huber's scale estimator under normality
  fct.Huber <- function(r) # central part
  {r^2*exp(-r^2/2)/(sqrt(2*pi))}    
  beta1= integrate(fct.Huber,lower=-k,upper=k)$value
  beta2=k^2*(pnorm(-k)+1-pnorm(k))
  beta=beta1+beta2
  return(beta)
}

invlogit<-function(u){return(1/(1+exp(-u)))}


#############################################
###########  linear regression  #############
#############################################

clipped_linear <- function(x,y,k=1.345,fisher_beta=0.7101645,scale=T,private=T,mu=1,maxiter=100,eps=1e-5,beta0=rep(0,dim(x)[2]),s0=1,stepsize=NULL,s.min=0.01,stopping=0,c_grad)
{
  #y is a vector of length n
  #x is an n by p matrix
  n=length(x[,1])
  p=length(x[1,])
  
  if(k!=1.345)            fisher_beta = Fisher.constant(k)
  
  noise=0 #gets updated further down if private=T
  conv=F
  iter=0
  grad=1

  
  r=(y-as.vector(x%*%beta0))/s0
  psi.vec<-psiHuber(r,k)
  
  rownorms<-apply(x,1,norm,type="2") #this does not need to get updated at each iteration

  
  # Location estimation with known scale
  if(scale==F)
  {
    eta<-ifelse(is.null(stepsize),(1/2)-0.001,stepsize)
    if(private==T) {
      noise=2*c_grad/(n*(mu/sqrt(maxiter))) #when scale is known, global sensitivity is 2*c_grad/n
      eps=max(eps,noise)
    }
    
    s=s0
    
    weightvec<-1/pmax(1,rownorms*abs(psi.vec)/c_grad) #this WILL get updated at each iteration
    
    # this performs gradient descent to estimate beta (if private=T, the gradient descent is noisy)
    while(iter<maxiter & grad>stopping*eps){
      iter=iter+1 
      
      beta=beta0+eta*(colMeans(psi.vec*weightvec*x)+noise*rnorm(p))
      beta0=beta
      
      r=(y-as.vector(x%*%beta0))/s0
      psi.vec<-psiHuber(r,k)
      
      #update the weights (with clipping, they must be recalculated at every iteration)
      weightvec<-1/pmax(1,rownorms*abs(psi.vec)/c_grad)
      
      grad=sqrt(sum((colMeans(psi.vec*weightvec*x)^2))) #this is the actual gradient evaluated at current value of beta0, NOT a noisy copy
    }
    beta=beta0
  }
  
  
  # Joint location and scale estimation
  if(scale==T)
  {
    theta0=c(beta0,s0)

    GS=2*c_grad  # global sensitivity in L2 norm
    eta<-ifelse(is.null(stepsize),(1/sqrt(1+k^2))-0.001,stepsize)   # step size 
    
    if(private==T)
    {
      noise=GS/(n*(mu/sqrt(maxiter))) 
      eps=max(eps,noise)
    }
    
    theta_norms<-sqrt((rownorms*abs(psi.vec))^2+((psiHuber(r,k)^2-fisher_beta)/2)^2)
    weightvec<-1/pmax(1,theta_norms/c_grad)
    sum.chi=mean(((psiHuber(r,k)^2)-fisher_beta)*weightvec)/2 ##divide by 2 so psi and chi come from same objective function
    
    while(iter<maxiter & grad > stopping*eps ) 
    {
      iter=iter+1 
      
      theta0=theta0+s0*eta*(c(colMeans(psi.vec*weightvec*x),sum.chi)+noise*rnorm(p+1)) #this is with eta*noise
      beta0=theta0[1:p]
      s0=theta0[(p+1)]

      r=(y-as.vector(x%*%beta0))/s0
      psi.vec<-psiHuber(r,k)
      
      theta_norms<-sqrt((rownorms*abs(psi.vec))^2+((psiHuber(r,k)^2-fisher_beta)/2)^2)
      weightvec<-1/pmax(1,theta_norms/c_grad)
      
      sum.chi=mean(((psiHuber(r,k)^2)-fisher_beta)*weightvec)/2
      
      grad=sqrt(sum((colMeans(psi.vec*weightvec*x)^2))+(sum.chi^2)) #non-noisy gradient evaluated at (beta0,sigma0)

      if(s0 < 0) grad=eps+1    # force the algorithm to continue if the current scale estimate is negative
    }
    
    beta=beta0
    s=s0
    
  }
  
  if( grad < eps) conv=T #keep in mind that this is the non-noisy gradient, so can use this to assess performance in
  #simulations, but NOT in real implementations (would violate privacy guarantee)
  
  out=NULL
  out$beta=beta
  out$s=s
  out$r=r
  out$iter=iter
  out$grad=grad #for implementation, need to change this to be the private version of gradient
  out$converge=1*conv #for implementation, need to change this to be based on the private version of gradient
  out$private=private
  out$noise=noise

  
  return(out)
}


#############################################
########### logistic regression #############
#############################################

clipped_logistic<-function(x,y,private=F,mu=1,maxiter,eps=1e-5,beta0=rep(0,dim(x)[2]),stopping=0,stepsize=NULL,c_grad){
  n=length(x[,1])
  p=length(x[1,])
  grad<-1
  iter<-0
  conv=F
  noise<-0
  
  yhat<-invlogit(as.vector(x%*%beta0))
  diffs<-y-yhat
  
  rownorms<-apply(x,1,norm,type="2")
  weightvec<-1/pmax(1,rownorms*abs(diffs)/c_grad)
  
  GS<-2*c_grad  #previously this was missing a factor of 2
  eta<-ifelse(is.null(stepsize),(4/(c_grad^2))-0.001,stepsize) # not intended to be a meaningful choice of step size--
                                                               # just keeping it similar to the non-clipping version
  
  if(private==T){noise<-GS/(n*(mu/sqrt(maxiter)))}
  eps=max(eps,noise)
  
  while(iter<maxiter & grad > stopping*eps){
    iter<-iter+1
    
    beta<-beta0+eta*(colMeans(diffs*weightvec*x)+noise*rnorm(p))
    beta0<-beta
    yhat<-invlogit(as.vector(x%*%beta0))
    diffs<-y-yhat
    
    weightvec<-1/pmax(1,rownorms*abs(diffs)/c_grad)
    grad=sqrt(sum(colMeans(diffs*weightvec*x)^2))
  }
  beta<-beta0
  
  
  if( grad < eps) conv=T
  
  out=NULL
  out$beta=beta
  out$iter=iter
  out$converge=1*conv
  out$private=private
  out$grad=grad
  
  return(out)
  
}


