###############
#version notes
#
# starts from v6
# which adds correction factor to sandwich variance, and uses 
# gaussian wigner matrix method ("analyze gauss") for BOTH M and Q
# 
# this one truncates the eigenvalues of Q if our noisy Q is not positive definite
# (post-processing, does not impact the privacy guarantee)


################################################
# Noisy gradient descent for linear regression with Huber-style loss function
# 


####
# auxiliary functions
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



weightfn<-function(x,max.norm=sqrt(2)){
  return(min(1,(max.norm/norm(x,type="2"))^2))
}



#################
# Main function

NGD.Huber <- function(x,y,k=1.345,fisher_beta=0.7101645,scale=T,private=T,mu=1,maxiter=100,eps=1e-5,beta0=rep(0,dim(x)[2]),s0=1,stepsize=NULL,s.min=0.01,stopping=0,mnorm=sqrt(2))
{
    #y is a vector of length n
    #x is now an n by p matrix
    n=length(x[,1])
    p=length(x[1,])
   
  if(k!=1.345)            fisher_beta = Fisher.constant(k)
    
  noise=0
  middle_noise=0 #gets updated further down if private=T
  outer_noise=0
  conv=F
  iter=0
  grad=1
  truncation<-0
  

  grad_traj<-rep(NA,maxiter+1)
  sigma_vec<-s0

  r=(y-as.vector(x%*%beta0))/s0
  psi.vec<-psiHuber(r,k)
  
  weightvec<-apply(x,1,weightfn,max.norm=mnorm)
  
  sum.chi=mean(((psiHuber(r,k)^2)-fisher_beta)*weightvec)/2 ##divide by 2 so psi and chi come from same objective function
  
  
  
  
# Location estimation with known scale
  if(scale==F)
    {
    eta<-ifelse(is.null(stepsize),(1/(2*k))-0.001,stepsize) # step size
    if(private==T) {
      noise=2*mnorm*k/(n*(mu/sqrt(maxiter+2)))
      eps=max(eps,noise)
    }
    
    s=s0
    
    grad_traj[1]<-sqrt(sum((colMeans(psi.vec*weightvec*x)^2))) #tracks the evolution of non-noisy gradient (in L2 norm)
    
    # this performs gradient descent to estimate beta (if private=T, the gradient descent is noisy)
  while(iter<maxiter & grad > stopping*eps)
    {
    iter=iter+1 
    
    beta=beta0+eta*(colMeans(psi.vec*weightvec*x)+noise*rnorm(p))
    beta0=beta

    r=(y-as.vector(x%*%beta0))/s0
    psi.vec<-psiHuber(r,k)
    
    
    grad=sqrt(sum((colMeans(psi.vec*weightvec*x)^2))) #this is the actual gradient evaluated at current value of beta0, NOT a noisy copy
    grad_traj[iter+1]<-grad
    
  }
  

  
    ## perform inference based on sandwich estimator
    
    ## compute estimate of "M" matrix and add noise to make private
    ## using analyze gauss (gaussian wigner matrix) mechanism for M privacy
    outer_term<-matrix(0,nrow=p,ncol=p)

    for(i in 1:n){
      outer_term<-outer_term+(abs(r[i])<k)*weightvec[i]*(x[i,]%o%x[i,])
    }
    
    outer_term<-outer_term/n #outer term as a matrix, we want the average
  

    if(private==T) {outer_noise<-(mnorm^2)/(n*(mu/sqrt(maxiter+2)))}
    
    outer_noisevec<-outer_noise*rnorm(p*(p-1)/2)
    
    outer_noisematrix<-matrix(0,nrow=p,ncol=p)
    outer_noisematrix[upper.tri(outer_noisematrix,diag=FALSE)]<-outer_noisevec
    
    #reflect them across the diagonal to get a symmetric matrix,
    #and draw p-many more random normals to put on the diagonal itself
    outer_noisematrix<-outer_noisematrix+t(outer_noisematrix)+diag(x=outer_noise*rnorm(p),nrow=p)
    
    outer_term<-outer_term+outer_noisematrix
    
    outer_term<-outer_term/s0 #finally, divide through by s0 (private estimate of it)

    
    ## next compute estimate of "Q" matrix and add noise  to make private
    
    
    #draw a vector of (appropriately scaled) random normal variables
    if(private==T) {middle_noise<-(mnorm^2)*(k^2)/(n*(mu/sqrt(maxiter+2)))}
    
    noisevec<-middle_noise*rnorm(p*(p-1)/2)
    
    #arrange them into the upper triangle of a matrix (leaving diagonal blank for now)
    noisematrix<-matrix(0,nrow=p,ncol=p)
    noisematrix[upper.tri(noisematrix,diag=FALSE)]<-noisevec
    
    #reflect them across the diagonal to get a symmetric matrix,
    #and draw p-many more random normals to put on the diagonal itself
    noisematrix<-noisematrix+t(noisematrix)+diag(x=middle_noise*rnorm(p),nrow=p)
    
    middle_term<-matrix(0,nrow=p,ncol=p)
    
    for(i in 1:n){
      middle_term<-middle_term+(psiHuber(r[i],k)^2)*((weightvec[i])^2)*(x[i,]%o%x[i,])
    }
    
    middle_term<-middle_term/n #we want an average
    
    middle_term<-middle_term+noisematrix #add the noise matrix to make private
    
    ## by definition Q should be positive definite; if adding noise causes
    ## any eigenvalues to become negative, truncate those eigenvalues
    ## to a small positive number
    decomp<-eigen(middle_term)
    
    if(min(decomp$values)<0){
      truncation<-1
      lambda<-decomp$values
      lambda[lambda<0]<-min(1/n,0.0001)
      
      middle_term<-(decomp$vectors)%*%diag(lambda)%*%solve(decomp$vectors)
    }
    
    sandwich<-solve(outer_term)%*%middle_term%*%t(solve(outer_term))
    
    corrected_variances<-(sandwich/n)+(noise^2)*diag(p)*(eta^2) #values on the diagonal of this matrix are corrected variances for the components of the beta vector
  
  }
  
  
  # Joint location and scale estimation
  if(scale==T)
    {
    theta0=c(beta0,s0)
     
    GS=sqrt(4*(mnorm^2)*k^2+(k^4)/4)
    eta<-ifelse(is.null(stepsize),(1/sqrt(1+k^2))-0.001,stepsize)             # step size 

    
    if(private==T)
      {
        noise=GS/(n*(mu/sqrt(maxiter+2))) 
        eps=max(eps,noise/2)
      }
    
    grad_traj[1]<-sqrt(sum((colMeans(psi.vec*weightvec*x)^2))+(sum.chi^2))
    
    while(iter<maxiter & grad > stopping*eps ) 
    {
      iter=iter+1 
      
      theta0=theta0+s0*eta*(c(colMeans(psi.vec*weightvec*x),sum.chi)+noise*rnorm(p+1)) #this is with eta*noise
      
      beta0=theta0[1:p]
      
      
      s0=theta0[(p+1)]
      sigma_vec<-c(sigma_vec,s0) #tracks the evolution of the scale estimate over iterations
      
      r=(y-as.vector(x%*%beta0))/s0
      
      psi.vec<-psiHuber(r,k)
      

      sum.chi=mean(((psiHuber(r,k)^2)-fisher_beta)*weightvec)/2
      
      grad=sqrt(sum((colMeans(psi.vec*weightvec*x)^2))+(sum.chi^2)) #non-noisy gradient evaluated at (beta0,sigma0)

      grad_traj[iter+1]<-grad #tracks the evolution of non-noisy gradient (in L2 norm)
      
      if(s0 < 0) grad=eps+1    # force the algorithm to continue if the current scale estimate is negative
      

    }
    
    
    beta=beta0
    s=s0
    
    
    ## perform inference based on sandwich estimator
    
    ## compute estimate of "M" matrix and add noise to make private
    ## using analyze gauss (gaussian wigner matrix) mechanism for M privacy
    outer_term<-matrix(0,nrow=p,ncol=p)
    
    for(i in 1:n){
      outer_term<-outer_term+(abs(r[i])<k)*weightvec[i]*(x[i,]%o%x[i,])
    }
    
    outer_term<-outer_term/n #outer term as a matrix, we want the average
    
    
    if(private==T) {outer_noise<-(mnorm^2)/(n*(mu/sqrt(maxiter+2)))}
    
    outer_noisevec<-outer_noise*rnorm(p*(p-1)/2)
    
    outer_noisematrix<-matrix(0,nrow=p,ncol=p)
    outer_noisematrix[upper.tri(outer_noisematrix,diag=FALSE)]<-outer_noisevec
    
    #reflect them across the diagonal to get a symmetric matrix,
    #and draw p-many more random normals to put on the diagonal itself
    outer_noisematrix<-outer_noisematrix+t(outer_noisematrix)+diag(x=outer_noise*rnorm(p),nrow=p)
    
    outer_term<-outer_term+outer_noisematrix
    

    outer_term<-outer_term/s0 #finally, divide through by s0 (private estimate of it)
    
    ######## using analyze gauss (gaussian wigner matrix) for Q
    ####################
    #####################
    
    
    #draw a vector of (appropriately scaled) random normals
    if(private==T) {middle_noise<-(mnorm^2)*(k^2)/(n*(mu/sqrt(maxiter+2)))}
    
    noisevec<-middle_noise*rnorm(p*(p-1)/2)
    
    #arrange them into the upper triangle of a matrix (leaving diagonal blank for now)
    noisematrix<-matrix(0,nrow=p,ncol=p)
    noisematrix[upper.tri(noisematrix,diag=FALSE)]<-noisevec
    
    #reflect them across the diagonal to get a symmetric matrix,
    #and draw p-many more random normals to put on the diagonal itself
    noisematrix<-noisematrix+t(noisematrix)+diag(x=middle_noise*rnorm(p),nrow=p)
    
    middle_term<-matrix(0,nrow=p,ncol=p)
    
    for(i in 1:n){
      middle_term<-middle_term+(psiHuber(r[i],k)^2)*((weightvec[i])^2)*(x[i,]%o%x[i,])
    }
    
    middle_term<-middle_term/n #we want an average
    
    middle_term<-middle_term+noisematrix #add the noise matrix to make private
    
    #truncate any negative eigenvalues of the Q matrix, deterministically (postprocessing, does not degrade privacy guarantee)
    decomp<-eigen(middle_term)
    
    if(min(decomp$values)<0){
      truncation<-1
      lambda<-decomp$values
      lambda[lambda<0]<-min(1/n,0.0001)
      
      middle_term<-(decomp$vectors)%*%diag(lambda)%*%solve(decomp$vectors)
    }
    
    
    sandwich<-solve(outer_term)%*%middle_term%*%t(solve(outer_term))
  
    corrected_variances<-(sandwich/n)+(noise^2)*diag(p)*(eta^2)
    
    
  }
  
  if( grad < eps) conv=T #keep in mind that this is the non-noisy gradient, so can use this to assess performance in
                          #simulations, but NOT in real implementations (would violate privacy guarantee)
    
  out=NULL
  out$beta=beta
  out$s=s
  out$r=r
  out$sigma_vec=sigma_vec
  out$iter=iter
  out$grad=grad #for implementation, need to change this to be the private version of gradient
  out$converge=1*conv #for implementation, need to change this to be based on the private version of gradient
  out$private=private
  out$noise=noise
  out$middle=middle_term
  out$outer=outer_term
  out$sandwich=sandwich
  out$variances=corrected_variances
  out$truncation=truncation
  out$gradtraj=grad_traj
  

  return(out)
}
