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

NoisyNewton <- function(x,y,k=1.345,fisher_beta=0.7101645,scale=T,private=T,mu=1,maxiter=100,eps=1e-10,beta0=rep(0,length(x[1,])),s0=1,s.min=0.01,stepsize=1,stopping=0,suppress.inference=FALSE,mnorm=sqrt(2)){
  
  #y is a vector of length n
  #x is an n by p matrix
  n=length(x[,1])
  p=length(x[1,])
  
  if(k!=1.345)            fisher_beta = Fisher.constant(k)
  
  noise=grad_noise=other_noise=middle_noise=outer_noise=hessian_noise=0 #gets updated further down if private=T
  conv_np=conv_priv=F
  iter=0
  grad=1
  middle_term<-NA
  outer_term<-NA
  sandwich<-NA
  var_correction<-NA
  stop_flag<-0
  eig_min<-NA
  
  np_grad_traj=priv_grad_traj=rep(NA,maxiter+1)
  
  weightvec<-apply(x,1,weightfn,max.norm=mnorm)
  
  r=(y-as.vector(x%*%beta0))/s0
  psi.vec<-psiHuber(r,k)
  psi2<-ifelse(abs(psi.vec)<k,psi.vec,0)
  sum.chi=mean(((psiHuber(r,k)^2)-fisher_beta)*weightvec)/2 ##divide by 2 so psi and chi come from same objective function
  
  
  sum.chiprime=mean((r^2)*(abs(r)<k)*weightvec)
  
  


  
  
  # Location estimation with known scale
  if(scale==F){
    
    if(private==T){
      grad_noise=2*mnorm*k/(n*(mu/sqrt(2*maxiter+2))) 
      #we need to do maxiter steps for the gradient AND maxiter steps for the Hessian, plus 2 for M and Q matrices
      #so sqrt(2*maxiter+2) in each place we add noise
      hessian_noise<-(mnorm^2)/(n*(mu/sqrt(2*maxiter+2))) #gets calculated at every step
      outer_noise<-(mnorm^2)/(n*(mu/sqrt(2*maxiter+2))) #calculated only once
      middle_noise<-(mnorm^2)*(k^2)/(n*(mu/sqrt(2*maxiter+2))) #calculated only once
      
      eps=max(eps,grad_noise/2)
    }
    
    s=s0
    
    np_grad_traj[1]<-sqrt(sum((colMeans(psi.vec*weightvec*x)^2))) #tracks the evolution of non-noisy gradient (in L2 norm)
    
    noisy_grad<-colMeans(psi.vec*weightvec*x)+grad_noise*rnorm(p)
    priv_grad_traj[1]<-sqrt(sum(noisy_grad^2))
    
    #initialize the hessian before starting newton iterations
    beta_hessian<-matrix(0,nrow=p,ncol=p)
    
    for(i in 1:n){
      beta_hessian<-beta_hessian+(abs(r[i])<k)*weightvec[i]*(x[i,]%o%x[i,])
    }
    
    beta_hessian<-beta_hessian/n
    
    
    ##sample the noise for the hessian (which will just be 0's if private=F)
    hessian_noisevec<-hessian_noise*rnorm(p*(p-1)/2)
    hessian_noisematrix<-matrix(0,nrow=p,ncol=p)
    hessian_noisematrix[upper.tri(hessian_noisematrix,diag=FALSE)]<-hessian_noisevec
    ##reflect them across the diagonal to get a symmetric matrix,
    ##and draw p-many more random normals to put on the diagonal itself
    hessian_noisematrix<-hessian_noisematrix+t(hessian_noisematrix)+diag(x=hessian_noise*rnorm(p),nrow=p)
    
    noisy_hessian<-beta_hessian+hessian_noisematrix
   

    
    #check whether stopping conditions are met
    if(stopping=="private" & priv_grad_traj[1]<=eps){
      stop_flag<-1
    }else if(stopping=="non-private" & np_grad_traj[1]<=eps){
      stop_flag<-1
    }
    # this performs newton raphson to estimate beta
    while(iter < maxiter & stop_flag < 1 & min(abs(eigen(noisy_hessian)$values))>1e-15 ){ #stop if hessian is computationally singular
    
      iter=iter+1
    
      beta=beta0+stepsize*solve((noisy_hessian)/s0)%*%(noisy_grad)
      beta0=beta
      
      r=(y-as.vector(x%*%beta0))/s0
      psi.vec<-psiHuber(r,k)

      #reset the hessian

      beta_hessian<-matrix(0,nrow=p,ncol=p)
      for(i in 1:n){
        beta_hessian<-beta_hessian+(abs(r[i])<k)*weightvec[i]*(x[i,]%o%x[i,])
      }
      beta_hessian<-beta_hessian/n
      
      ##sample the noise for the hessian
      hessian_noisevec<-hessian_noise*rnorm(p*(p-1)/2)
      hessian_noisematrix<-matrix(0,nrow=p,ncol=p)
      hessian_noisematrix[upper.tri(hessian_noisematrix,diag=FALSE)]<-hessian_noisevec
      hessian_noisematrix<-hessian_noisematrix+t(hessian_noisematrix)+diag(x=hessian_noise*rnorm(p),nrow=p)
      
      noisy_hessian<-beta_hessian+hessian_noisematrix
      
      grad=sqrt(sum((colMeans(psi.vec*weightvec*x)^2))) #this is the actual (norm of the) gradient evaluated at current value of beta0, NOT a noisy copy
      np_grad_traj[iter+1]<-grad
      noisy_grad<-colMeans(psi.vec*weightvec*x)+grad_noise*rnorm(p)
      priv_grad_traj[iter+1]<-sqrt(sum(noisy_grad^2))
      
      #check whether stopping conditions are met
      if(stopping=="private" & priv_grad_traj[iter+1]<=eps){
          stop_flag<-1
      }else if(stopping=="non-private" & np_grad_traj[iter+1]<=eps){
          stop_flag<-1
      }
      
    }
    
    beta=beta0
    
  }
  
  
  # Joint location and scale estimation
  
  if(scale==T)
  {
    theta0=c(beta0,s0)
    
    GS_sigmagrad<-min((k^2)/2,(k^2-fisher_beta))
    GS_grad<-sqrt(4*(mnorm^2)*k^2+GS_sigmagrad^2)
    if(private==T){
      grad_noise<-GS_grad/(n*(mu/sqrt(3*maxiter+2)))
      hessian_noise<-2/(n*(mu/sqrt(3*maxiter+2)))
      other_noise<-sqrt(4*(mnorm^2)*(k^2)+(k^4))/(n*(mu/sqrt(3*maxiter+2)))
      outer_noise<-(mnorm^2)/(n*(mu/sqrt(3*maxiter+2)))
      middle_noise<-(mnorm^2)*(k^2)/(n*(mu/sqrt(3*maxiter+2)))
      eps=max(eps,grad_noise/2)
    }
    
    grad<-sqrt(sum((colMeans(psi.vec*weightvec*x)^2))+(sum.chi^2))
    np_grad_traj[1]<-grad
    
    noisy_grad<-c(colMeans(psi.vec*weightvec*x),sum.chi)+grad_noise*rnorm(p+1)
    priv_grad_traj[1]<-sqrt(sum(noisy_grad^2))
    
    
    #initialize the hessian before starting newton iterations
    beta_hessian<-matrix(0,nrow=p,ncol=p)
    
    for(i in 1:n){
      beta_hessian<-beta_hessian+(abs(r[i])<k)*weightvec[i]*(x[i,]%o%x[i,])
    }
    
    beta_hessian<-beta_hessian/n
    
    
    ##sample the noise for the hessian (which will just be 0's if private=F)
    hessian_noisevec<-hessian_noise*rnorm(p*(p-1)/2)
    hessian_noisematrix<-matrix(0,nrow=p,ncol=p)
    hessian_noisematrix[upper.tri(hessian_noisematrix,diag=FALSE)]<-hessian_noisevec
    ##reflect them across the diagonal to get a symmetric matrix,
    ##and draw p-many more random normals to put on the diagonal itself
    hessian_noisematrix<-hessian_noisematrix+t(hessian_noisematrix)+diag(x=hessian_noise*rnorm(p),nrow=p)
    
    noisy_hessian<-matrix(nrow=p+1,ncol=p+1)
    noisy_hessian[1:p,1:p]<-beta_hessian+hessian_noisematrix
    
    mixed_partial<-colMeans(psi2*weightvec*x) #the d(beta)d(sigma) elements of the full hessian
    lastrow<-c(mixed_partial,sum.chiprime)+other_noise*rnorm(p+1)
    
    noisy_hessian[,(p+1)]<-lastrow
    noisy_hessian[(p+1),]<-lastrow
    
    true_hessian<-matrix(nrow=p+1,ncol=p+1)
    true_hessian[1:p,1:p]<-beta_hessian
    true_hessian[,(p+1)]<-c(mixed_partial,sum.chiprime)
    true_hessian[(p+1),]<-c(mixed_partial,sum.chiprime)
    
    #check whether stopping conditions are met
    if(stopping=="private" & priv_grad_traj[1]<=eps){
        stop_flag<-1
    }else if(stopping=="non-private" & np_grad_traj[1]<=eps){
        stop_flag<-1
    }
    
    eig_min<-min(eigen(true_hessian)$value)
    
    # this performs newton raphson to estimate beta
    while(iter < maxiter & stop_flag < 1){ 
      
      iter=iter+1
      
      eig_min<-min(eigen(true_hessian)$value)
      
      theta=theta0+stepsize*solve(noisy_hessian/s0)%*%noisy_grad
      theta0=theta
      
      beta0=theta0[1:p]
      s0=theta0[(p+1)]
      
      r=(y-as.vector(x%*%beta0))/s0
      psi.vec<-psiHuber(r,k)
      psi2<-ifelse(abs(psi.vec)<k,psi.vec,0)
      
      sum.chi=mean(((psiHuber(r,k)^2)-fisher_beta)*weightvec)/2
      sum.chiprime=mean((r^2)*(abs(r)<k)*weightvec)
      
      grad<-sqrt(sum((colMeans(psi.vec*weightvec*x)^2))+(sum.chi^2))
      np_grad_traj[iter+1]<-grad
      
      noisy_grad<-c(colMeans(psi.vec*weightvec*x),sum.chi)+grad_noise*rnorm(p+1)
      priv_grad_traj[iter+1]<-sqrt(sum(noisy_grad^2))
      
      #reset the hessian
      beta_hessian<-matrix(0,nrow=p,ncol=p)
      for(i in 1:n){
        beta_hessian<-beta_hessian+(abs(r[i])<k)*weightvec[i]*(x[i,]%o%x[i,])
      }
      beta_hessian<-beta_hessian/n
      
      ##sample the noise for the hessian
      hessian_noisevec<-hessian_noise*rnorm(p*(p-1)/2)
      hessian_noisematrix<-matrix(0,nrow=p,ncol=p)
      hessian_noisematrix[upper.tri(hessian_noisematrix,diag=FALSE)]<-hessian_noisevec
      hessian_noisematrix<-hessian_noisematrix+t(hessian_noisematrix)+diag(x=hessian_noise*rnorm(p),nrow=p)
      
      true_hessian<-matrix(nrow=p+1,ncol=p+1)
      true_hessian[1:p,1:p]<-beta_hessian
      
      noisy_hessian<-matrix(nrow=p+1,ncol=p+1)
      noisy_hessian[1:p,1:p]<-beta_hessian+hessian_noisematrix
      
      mixed_partial<-colMeans(psi2*weightvec*x) #the d(beta)d(sigma) elements of the full hessian
      lastrow<-c(mixed_partial,sum.chiprime)+other_noise*rnorm(p+1)
      
      true_hessian[,(p+1)]<-c(mixed_partial,sum.chiprime)
      true_hessian[(p+1),]<-c(mixed_partial,sum.chiprime)
      
      noisy_hessian[,(p+1)]<-lastrow
      noisy_hessian[(p+1),]<-lastrow
      

      
      
      #check whether stopping conditions are met
      if(stopping=="private" & priv_grad_traj[iter+1]<=eps){
          stop_flag<-1
      }else if(stopping=="non-private" & np_grad_traj[iter+1]<=eps){
          stop_flag<-1
      }
      
    }#end of while loop for iterations 
    beta=beta0
    s=s0
    
  }#end of "joint estimation" section
  
    if(suppress.inference==FALSE){
      ## perform inference based on sandwich estimator
      
      ## compute estimate of "M" matrix and add noise to make private
      ## using analyze gauss (gaussian wigner matrix) mechanism for M privacy
      outer_term<-matrix(0,nrow=p,ncol=p)
      
      for(i in 1:n){
        outer_term<-outer_term+(abs(r[i])<k)*weightvec[i]*(x[i,]%o%x[i,])
      }
      
      outer_term<-outer_term/n #outer term as a matrix, we want the average
      
      
      outer_noisevec<-outer_noise*rnorm(p*(p-1)/2)
      
      outer_noisematrix<-matrix(0,nrow=p,ncol=p)
      outer_noisematrix[upper.tri(outer_noisematrix,diag=FALSE)]<-outer_noisevec
      
      #reflect them across the diagonal to get a symmetric matrix,
      #and draw p-many more random normals to put on the diagonal itself
      outer_noisematrix<-outer_noisematrix+t(outer_noisematrix)+diag(x=outer_noise*rnorm(p),nrow=p)
      
      outer_term<-outer_term+outer_noisematrix
      
      outer_term<-outer_term/s0 #finally, divide through by s0 (private estimate of it)
      
      
      
      ## next compute estimate of "Q" matrix and add noise to make private
      
      
      #draw a vector of (appropriately scaled) random normal variables
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
      
      hessian_inverse<-solve(outer_term)
      
      sandwich<-hessian_inverse%*%middle_term%*%t(hessian_inverse)
      
      var_correction<-hessian_inverse%*%diag(grad_noise^2,nrow=p)%*%hessian_inverse
      
      corrected_variances<-(sandwich/n)+(var_correction)*stepsize^2
      

    }
    
    
  
  
  if(grad < eps) conv_np=T 
  if(sqrt(sum(noisy_grad^2)) < eps ) conv_priv=T 
  
  out=NULL
  out$beta=beta
  out$s=s
  out$r=r
  out$iter=iter
  out$grad=grad
  out$nonpriv_converge=1*conv_np #convergence assessment based on non-private version of gradient
  out$priv_converge=1*conv_priv #convergence assessment based on private version of gradient  out$private=private
  out$eig_min=eig_min
  out$noise=grad_noise
  out$hess<-beta_hessian
  out$priv_gradtraj=priv_grad_traj
  out$nonpriv_gradtraj=np_grad_traj
  if(suppress.inference==FALSE){
  out$middle=middle_term
  out$outer=outer_term
  out$sandwich=sandwich
  out$variances=corrected_variances
  }
  
  return(out)
}

