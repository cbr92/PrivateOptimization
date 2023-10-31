source("src/extra/Newton_linear_regression_extra.R")

n<-1000
beta_gen=c(1,1,1,1)
s0=2
ns<-500

stepsizes<-c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
diverged<-rep(NA,length(stepsizes))
converged<-rep(NA,length(stepsizes))

set.seed(85330)

for(j in 1:10){
  niter<-rep(NA,ns)
  conv<-rep(NA,ns)
  eigs<-rep(NA,ns)
  
  for(i in 1:ns){
    x=matrix(rnorm(n*4,mean=0,sd=2),nrow=n)
    x[,1]<-rep(1,n)
    y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
    out=NoisyNewton(x=x,y=y,private=T,mu=2,scale=T,stepsize=stepsizes[j],maxiter=20,suppress.inference=TRUE,stopping="private",mnorm=sqrt(2)) #vary the step size
    
    niter[i]<-out$iter
    conv[i]=out$priv_converge
    eigs[i]<-out$eig_min
  }
  
  diverged[j]<-mean(eigs<1e-14)
  converged[j]<-mean(conv)
}


plot(1:length(stepsizes),diverged,xlab="Step Size",ylab="% Trials Diverging",xaxt='n',pch=16, main=NULL)
axis(1, at=1:10, labels=stepsizes)

