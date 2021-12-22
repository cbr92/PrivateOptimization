# contains differentially private implementations of 
# gradient descent and newton's method for 
# logistic regression models
#
#############################################



#############################################
########### auxiliary functions #############
#############################################
invlogit<-function(u){return(1/(1+exp(-u)))}

hess<-function(v,xmat){
  u<-exp(-as.vector(xmat%*%v))
  return(u/((1+u)^2))
}

weightfn<-function(x,max.norm=sqrt(2)){
  return(min(1,(max.norm/norm(x,type="2"))^2))
}

#############################################
####### gradient descent function ###########
#############################################

NGD_logistic<-function(x,y,private=F,mu=1,maxiter,eps=1e-5,beta0=rep(0,dim(x)[2]),stopping=0,stepsize=NULL,mnorm=sqrt(2),suppress.inference=FALSE){
  n=length(x[,1])
  p=length(x[1,])
  
  iter<-0
  truncation<-0
  conv=F
  noise=outer_noise=middle_noise=0
  
  yhat<-invlogit(as.vector(x%*%beta0))
  diffs<-y-yhat
  weightvec<-apply(x,1,weightfn,max.norm=mnorm)
  
  GS<-2*mnorm
  eta<-ifelse(is.null(stepsize),(4/(mnorm^2))-0.001,stepsize) #step size
  if(private==T){
    noise<-GS/(n*(mu/sqrt(maxiter+2))) #the +2 is needed to get inference (M and Q) within the total mu privacy budget
    outer_noise<-((mnorm^2)/4)/(n*(mu/sqrt(maxiter+2)))
    middle_noise<-(mnorm^2)/(n*(mu/sqrt(maxiter+2)))
    eps=max(eps,noise)
  }
  
  noisy_grad<-colMeans(diffs*weightvec*x)+noise*rnorm(p)
  
  while(iter < maxiter & sqrt(sum(noisy_grad^2)) > stopping*eps){
    iter<-iter+1
    old_grad<-noisy_grad
    
    beta<-beta0+eta*noisy_grad
    beta0<-beta
    yhat<-invlogit(as.vector(x%*%beta0))
    diffs<-y-yhat
    
    noisy_grad<-colMeans(diffs*weightvec*x)+noise*rnorm(p)
  }
  beta<-beta0
  
  if(suppress.inference==FALSE){
    ## perform inference based on sandwich estimator
    
    ## compute estimate of "M" matrix and add noise to make private
    ## using analyze gauss (gaussian wigner matrix) mechanism for M privacy
    outer_term<-matrix(0,nrow=p,ncol=p)
    
    for(i in 1:n){
      outer_term<-outer_term+(yhat[i]*(1-yhat[i]))*weightvec[i]*(x[i,]%o%x[i,])
    }
    
    outer_term<-outer_term/n #outer term as a matrix, we want the average
        
    outer_noisevec<-outer_noise*rnorm(p*(p-1)/2)
    
    outer_noisematrix<-matrix(0,nrow=p,ncol=p)
    outer_noisematrix[upper.tri(outer_noisematrix,diag=FALSE)]<-outer_noisevec
    
    #reflect them across the diagonal to get a symmetric matrix,
    #and draw p-many more random normals to put on the diagonal itself
    outer_noisematrix<-outer_noisematrix+t(outer_noisematrix)+diag(x=outer_noise*rnorm(p),nrow=p)
    
    outer_term<-outer_term+outer_noisematrix
    
    
    ## next compute estimate of "Q" matrix and add noise  to make private
    
    
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
      middle_term<-middle_term+(diffs[i]^2)*((weightvec[i])^2)*(x[i,]%o%x[i,])
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
  final_grad<-ifelse(iter<maxiter,noisy_grad,old_grad)
  if( sqrt(sum(final_grad^2)) < eps) conv=T
  
  out=NULL
  out$beta=beta
  out$iter=iter
  out$converge=1*conv
  out$private=private
  out$grad=final_grad
  if(suppress.inference==FALSE){
    out$Q=middle_term
    out$M=outer_term
    out$sandwich=sandwich
    out$variances=corrected_variances
    out$truncation=truncation
  }
  
  return(out)
  
}



#############################################
####### newton's method function ############
#############################################



Newton_logistic<-function(x,y,private=F,mu=1,maxiter,eps=1e-5,beta0=rep(0,dim(x)[2]),stopping=0,stepsize=NULL,mnorm=sqrt(2),suppress.inference=FALSE){
  n=length(x[,1])
  p=length(x[1,])
  
  iter<-0
  truncation<-0
  conv=F
  noise=hessian_noise=outer_noise=middle_noise=0
  
  yhat<-invlogit(as.vector(x%*%beta0))
  diffs<-y-yhat
  hessian_coefs<-yhat*(1-yhat)
  weightvec<-apply(x,1,weightfn)
  
  eta<-ifelse(is.null(stepsize),1,stepsize) #perform pure newton unless otherwise specified
 
  if(private==T){
    noise<-2*sqrt(2)/(n*(mu/sqrt(2*maxiter+2))) #noise for the gradient
    hessian_noise<-((mnorm^2)/4)/(n*(mu/sqrt(2*maxiter+2))) #calculated at every iteration
    outer_noise<-((mnorm^2)/4)/(n*(mu/sqrt(2*maxiter+2))) #calculated only once
    middle_noise<-(mnorm^2)/(n*(mu/sqrt(2*maxiter+2)))  #calculated only once
    eps=max(eps,noise)
  }
  
  noisy_grad<-colMeans(diffs*weightvec*x)+noise*rnorm(p)
  
  beta_hessian<-matrix(0,nrow=p,ncol=p)
  for(i in 1:n){beta_hessian<-beta_hessian+weightvec[i]*hessian_coefs[i]*(x[i,]%o%x[i,])}
  beta_hessian<-beta_hessian/n 
  
  ##sample the noise for the hessian (which will just be 0's if private=F)
  hessian_noisevec<-hessian_noise*rnorm(p*(p-1)/2)
  hessian_noisematrix<-matrix(0,nrow=p,ncol=p)
  hessian_noisematrix[upper.tri(hessian_noisematrix,diag=FALSE)]<-hessian_noisevec
  ##reflect them across the diagonal to get a symmetric matrix,
  ##and draw p-many more random normals to put on the diagonal itself
  hessian_noisematrix<-hessian_noisematrix+t(hessian_noisematrix)+diag(x=hessian_noise*rnorm(p),nrow=p)
    
  noisy_hessian<-beta_hessian+hessian_noisematrix
  
  while(iter < maxiter & sqrt(sum(noisy_grad^2)) > stopping*eps & min(abs(eigen(noisy_hessian)$values))>1e-15){ #stop if hessian is computationally singular
    iter<-iter+1
    
    old_grad<-noisy_grad 
    
    beta<-beta0+eta*solve(noisy_hessian)%*%(noisy_grad)
    beta0<-beta
    
    yhat<-invlogit(as.vector(x%*%beta0))
    diffs<-y-yhat
    hessian_coefs<-yhat*(1-yhat)
    
    
    #reset the hessian
    beta_hessian<-matrix(0,nrow=p,ncol=p)
    for(i in 1:n){beta_hessian<-beta_hessian+weightvec[i]*hessian_coefs[i]*(x[i,]%o%x[i,])}
    beta_hessian<-beta_hessian/n
    
    ##sample the noise for the hessian (which will just be 0's if private=F)
    hessian_noisevec<-hessian_noise*rnorm(p*(p-1)/2)
    hessian_noisematrix<-matrix(0,nrow=p,ncol=p)
    hessian_noisematrix[upper.tri(hessian_noisematrix,diag=FALSE)]<-hessian_noisevec
    ##reflect them across the diagonal to get a symmetric matrix,
    ##and draw p-many more random normals to put on the diagonal itself
    hessian_noisematrix<-hessian_noisematrix+t(hessian_noisematrix)+diag(x=hessian_noise*rnorm(p),nrow=p)
    
    noisy_hessian<-beta_hessian+hessian_noisematrix
    
    noisy_grad<-colMeans(diffs*weightvec*x)+noise*rnorm(p)
    
  }
  beta<-beta0
  
  if(suppress.inference==FALSE){
    ## perform inference based on sandwich estimator
    
    ## compute estimate of "M" matrix and add noise to make private
    ## using analyze gauss (gaussian wigner matrix) mechanism for M privacy
       
    outer_noisevec<-outer_noise*rnorm(p*(p-1)/2)
    
    outer_noisematrix<-matrix(0,nrow=p,ncol=p)
    outer_noisematrix[upper.tri(outer_noisematrix,diag=FALSE)]<-outer_noisevec
    
    #reflect them across the diagonal to get a symmetric matrix,
    #and draw p-many more random normals to put on the diagonal itself
    outer_noisematrix<-outer_noisematrix+t(outer_noisematrix)+diag(x=outer_noise*rnorm(p),nrow=p)
    
    outer_term<-beta_hessian+outer_noisematrix
    
    
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
      middle_term<-middle_term+(diffs[i]^2)*((weightvec[i])^2)*(x[i,]%o%x[i,])
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
    corrected_variances<-(sandwich/n)+(hessian_inverse%*%diag(noise^2,nrow=p)%*%hessian_inverse)*(eta^2)
    
  }
  final_grad<-ifelse(iter<maxiter,noisy_grad,old_grad)
  if( sqrt(sum(final_grad^2)) < eps) conv=T
  
  out=NULL
  out$beta=beta
  out$iter=iter
  out$converge=1*conv
  out$private=private
  out$grad=final_grad
  if(suppress.inference==FALSE){
    out$Q=middle_term
    out$M=outer_term
    out$sandwich=sandwich
    out$variances=corrected_variances
    out$truncation=truncation
  }
  
  return(out)
  
}
