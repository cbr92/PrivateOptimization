# contains differentially private implementations of 
# gradient descent and fisher scoring for 
# logistic regression models

# note: currently this only covers estimation,
# but does not include inference functionality

#############################################



#############################################
########### auxiliary functions #############
#############################################
invlogit<-function(u){return(1/(1+exp(-u)))}

hess<-function(v,xmat){
  u<-exp(-as.vector(xmat%*%v))
  return(u/((1+u)^2))
}

weightfn<-function(x){
  return(min(1,2/(norm(x,type="2")^2)))
}

#############################################
####### gradient descent function ###########
#############################################

NGD_logistic<-function(x,y,private=F,mu=1,maxiter,eps=1e-5,beta0=rep(0,dim(x)[2]),stopping=0){
  n=length(x[,1])
  p=length(x[1,])
  grad<-1
  iter<-0
  conv=F
  noise<-0
  
  yhat<-invlogit(as.vector(x%*%beta0))
  diffs<-y-yhat
  weightvec<-apply(x,1,weightfn)
  
  GS<-sqrt(2)
  eta<-1 ### PLACEHOLDER-- think about what step size is appropriate
  if(private==T){noise<-GS/(n*(mu/sqrt(maxiter+2)))} #the +2 is needed to get inference (M and Q) within the total mu privacy budget
  eps=max(eps,noise/2)
  
  while(grad>stopping*eps & iter<maxiter){
    iter<-iter+1
    
    
    beta<-beta0+eta*(colMeans(diffs*weightvec*x)+noise*rnorm(p))
    beta0<-beta
    yhat<-invlogit(as.vector(x%*%beta0))
    diffs<-y-yhat
    
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



#############################################
####### newton's method function ############
####### (actually Fisher scoring) ###########
#############################################



Newton_logistic<-function(x,y,private=F,mu=1,maxiter,eps=1e-5,beta0=rep(0,dim(x)[2]),stopping=0){
  n=length(x[,1])
  p=length(x[1,])
  grad<-1
  iter<-0
  conv=F
  noise<-0
  hessian_noise<-0
  
  yhat<-invlogit(as.vector(x%*%beta0))
  diffs<-y-yhat
  hessian_coefs<-hess(beta0,x)
  weightvec<-apply(x,1,weightfn)
  
 
  if(private==T){
    noise<-sqrt(2)/(n*(mu/sqrt(2*maxiter)))
    hessian_noise<-(1/2)/(n*(mu/sqrt(2*maxiter)))
    }
  eps=max(eps,noise/2)
  
  beta_hessian<-matrix(0,nrow=p,ncol=p)
  for(i in 1:n){beta_hessian<-beta_hessian+weightvec[i]*hessian_coefs[i]*(x[i,]%o%x[i,])}
  beta_hessian<-beta_hessian/n 
  
  while(grad>stopping*eps & iter<maxiter & min(eigen(beta_hessian)$values)>1e-15){
    iter<-iter+1
    
    ##sample the noise for the hessian (which will just be 0's if private=F)
    hessian_noisevec<-hessian_noise*rnorm(p*(p-1)/2)
    hessian_noisematrix<-matrix(0,nrow=p,ncol=p)
    hessian_noisematrix[upper.tri(hessian_noisematrix,diag=FALSE)]<-hessian_noisevec
    ##reflect them across the diagonal to get a symmetric matrix,
    ##and draw p-many more random normals to put on the diagonal itself
    hessian_noisematrix<-hessian_noisematrix+t(hessian_noisematrix)+diag(x=hessian_noise*rnorm(p),nrow=p)
    
    noisy_hessian<-beta_hessian+hessian_noisematrix
    
    #### QUESTION: does truncating the eigenvalues of the noisy hessian help with the divergence problem?
    decomp<-eigen(noisy_hessian)
    if(min(decomp$values)<0){
      truncation<-1
      lambda<-decomp$values
      lambda[lambda<0]<-min(1/n,0.001)
      
      noisy_hessian<-(decomp$vectors)%*%diag(lambda)%*%solve(decomp$vectors)
    }
    
    beta<-beta0+solve(noisy_hessian)%*%(colMeans(diffs*weightvec*x)+noise*rnorm(p))
    beta0<-beta
    yhat<-invlogit(as.vector(x%*%beta0))
    diffs<-y-yhat
    hessian_coefs<-hess(beta0,x)
    
    
    #reset the hessian
    beta_hessian<-matrix(0,nrow=p,ncol=p)
    for(i in 1:n){beta_hessian<-beta_hessian+weightvec[i]*hessian_coefs[i]*(x[i,]%o%x[i,])}
    beta_hessian<-beta_hessian/n
    
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
