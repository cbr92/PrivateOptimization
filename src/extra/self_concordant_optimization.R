#############################################
########## CONTENTS #########################
#############################################

## 1. auxiliary functions
## 2. noisy gradient descent with self-concordant pseudo-huber loss
## 3. noisy newton's method with self-concordant pseudo-huber loss



#############################################
##########    1.    #########################
#############################################

# Here we use a self-concordant pseudo-huber loss function proposed in
# Ostrovskii and Bach (2021) (Eqn 21 in that paper)

## D. M. Ostrovskii and F. Bach. Finite-sample analysis of M-estimators using self-concordance.
## Electronic Journal of Statistics, 15(1):326-391, 2021.

phi.pseudo<-function(v,c=1){
  v<-v/c
  r<-sqrt(1+4*v^2)
  calc<-(r-1+log((r-1)/(2*v^2)))/2
  calc[v==0]<-0
  
  return((c^2)*calc)
}


# first derivative of loss function
phiprime.pseudo<-function(v,c=1){
  v<-v/c
  r<-sqrt(1+4*v^2)
  calc<-2*v/(r-1)-1/v
  calc[v==0]<-0
  calc[abs(v)<0.00002]<-0
  return(c*calc)
}

# second derivative of loss function
phi2prime.pseudo<-function(v,c=1){
  v<-v/c
  r<-sqrt(1+4*v^2)
  calc<-2/(r*(1-r))+1/(v^2)
  calc[v==0]<-1
  calc[abs(calc)>1]<-1
  return(calc)
  
}



weightfn<-function(x){
  return(min(1,2/(norm(x,type="2")^2)))
}


#############################################
##########    2.    #########################
#############################################


SCNGD <- function(x,y,tau=1,scale=T,private=T,mu=1,maxiter=100,eps=1e-5,beta0=rep(0,dim(x)[2]),s0=1,stopping=0)
{
  #y is a vector of length n
  #x is now an n by p matrix
  n=length(x[,1])
  p=length(x[1,])


  noise=0 #gets updated further down if private=T
  conv=F
  iter=0
  grad_traj<-rep(NA,maxiter+1)
  
  r=(y-as.vector(x%*%beta0))/s0
  psi.vec<-phiprime.pseudo(r,c=tau)
  
  
  weightvec<-apply(x,1,weightfn)
  
  
  # Location estimation with known scale
  if(scale==F)
  {
    eta=(1/(2*tau))-0.001 # step size
    if(private==T) {
      noise=2*sqrt(2)*tau/(n*(mu/sqrt(maxiter))) 
      eps=max(eps,noise)
    }
    
    s=s0
    
    grad_traj[1]<-sqrt(sum((colMeans(psi.vec*weightvec*x)^2)))
    grad<-grad_traj[1]
    
    # this performs gradient descent to estimate beta (if private=T, the gradient descent is noisy)
    while(iter<maxiter & grad > stopping*eps){
      iter=iter+1 
      
      beta=beta0+eta*(colMeans(psi.vec*weightvec*x)+noise*rnorm(4))
      beta0=beta
      
      r=(y-as.vector(x%*%beta0))/s0
      psi.vec<-phiprime.pseudo(r,c=tau)
      
      grad_traj[iter+1]<-sqrt(sum((colMeans(psi.vec*weightvec*x)^2)))
      grad<-grad_traj[iter+1]
    }
    beta=beta0
    
  }
  
  
  # Joint location and scale estimation not currently included

  
  if( grad < eps) conv=T

  out=NULL
  out$beta=beta
  out$s=s
  out$r=r
  out$iter=iter
  out$grad=grad
  out$converge=1*conv
  out$private=private
  out$noise=noise
  out$gradtraj=grad_traj
  
  
  
  return(out)
}



#############################################
##########    3.    #########################
#############################################


SCNewton <- function(x,y,tau=1,scale=T,private=T,mu=1,maxiter=100,max.s.iter=25,eps=1e-5,s.eps=1e-3,beta0=rep(0,4),s0=1,s.min=0.01,stepsize=1,eta=1,stopping=0){
  
  #y is a vector of length n
  #x is now an n by p matrix
  n=length(x[,1])
  p=length(x[1,])
  

  #"noise" values get updated later on if private=T
  noise=0
  hessian_noise=0
  conv=F
  iter=0
  grad_traj<-rep(NA,maxiter+1)
  
  
  
  r=(y-as.vector(x%*%beta0))/s0
  psi.vec<-phiprime.pseudo(r,c=tau)
 
  weightvec<-apply(x,1,weightfn)
  

  # Location estimation with known scale
  if(scale==F)
  {
    if(private==T) {
      
      noise=2*sqrt(2)*tau/(n*(mu/sqrt(2*maxiter)))
      hessian_noise<-2/(n*(mu/sqrt(2*maxiter)))
      eps=max(eps,noise)
    }
    
    s=s0
    
    grad_traj[1]<-sqrt(sum((colMeans(psi.vec*weightvec*x)^2)))
    grad<-grad_traj[1]
    
    #initialize the hessian before starting newton iterations
    beta_hessian<-matrix(0,nrow=p,ncol=p)
    
    for(i in 1:n){
      beta_hessian<-beta_hessian+phi2prime.pseudo(r[i],c=tau)*weightvec[i]*(x[i,]%o%x[i,])
    }
    
    beta_hessian<-beta_hessian/n
    
  
    # this performs newton raphson to estimate beta (if private=T, the newton raphson is noisy)
    
    while(iter < maxiter & grad > stopping*eps ){
    
      iter=iter+1
      
      ##sample the noise for the hessian (which will just be 0's if private=F)
      hessian_noisevec<-hessian_noise*rnorm(p*(p-1)/2)
      hessian_noisematrix<-matrix(0,nrow=p,ncol=p)
      hessian_noisematrix[upper.tri(hessian_noisematrix,diag=FALSE)]<-hessian_noisevec
      ##reflect them across the diagonal to get a symmetric matrix,
      ##and draw p-many more random normals to put on the diagonal itself
      hessian_noisematrix<-hessian_noisematrix+t(hessian_noisematrix)+diag(x=hessian_noise*rnorm(p),nrow=p)
      

      noisy_grad<-colMeans(psi.vec*weightvec*x)+noise*rnorm(p)
      beta=beta0+stepsize*solve((beta_hessian+hessian_noisematrix)/s0)%*%(noisy_grad)
      beta0=beta
      
      r=(y-as.vector(x%*%beta0))/s0
      psi.vec<-phiprime.pseudo(r,c=tau)
      
      
      #reset the hessian
      beta_hessian<-matrix(0,nrow=p,ncol=p)
      for(i in 1:n){
        beta_hessian<-beta_hessian+phi2prime.pseudo(r[i],c=tau)*weightvec[i]*(x[i,]%o%x[i,])
      }
      beta_hessian<-beta_hessian/n
      
      #checking whether the gradient is getting small, which is the goal of the optimization
      grad_traj[iter+1]<-sqrt(sum((colMeans(psi.vec*weightvec*x)^2)))
      grad<-grad_traj[iter+1]
      
    }
    
    beta=beta0
    
  }
  
  # Joint location and scale estimation not currently included
  
  
  if( grad < eps) conv=T
  
  out=NULL
  out$beta=beta
  out$s=s
  out$r=r
  out$iter=iter
  out$grad=grad
  out$converge=1*conv
  out$private=private
  out$noise=noise
  out$gradtraj=grad_traj
 
  
  return(out)
}


