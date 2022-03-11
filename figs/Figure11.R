source("src/extra/Newton_linear_regression_extra.R")

coverage_calculator<-function(v){
  return(mean(abs(v)<qnorm(0.975)))
}
library(MASS)

beta_gen=rep(1,4)
s0=2
p<-length(beta_gen)

M<-c(1,0.5,0,0.5,2,0,0,0,2)
M<-matrix(M,nrow=3,byrow=T)

means<-c(0,0,0)

samplesizes<-c(250,500,750,1000,1500,2000,4000)
niter<-c(15,16,17,17,18,18,19)
ns<-1000


raw_coverage_mat<-matrix(NA,nrow=length(samplesizes),ncol=p)
corrected_coverage_mat<-matrix(NA,nrow=length(samplesizes),ncol=p)
np_coverage_mat<-matrix(NA,nrow=length(samplesizes),ncol=p)


set.seed(303)
for(j in 1:length(samplesizes)){
  n<-samplesizes[j]
  iter.j<-niter[j]
  
  zscores_raw=zscores_corrected=zscores_np=matrix(nrow=ns,ncol=p)
  
  for(i in 1:ns){
  x<-mvrnorm(n=n,mu=means,Sigma=M)
  x<-cbind(rep(1,n),x)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  
  out_priv<-NoisyNewton(x=x,y=y,private=T,mu=1,scale=F,s0=2,beta0=rep(0,4),stepsize=0.50,maxiter=iter.j,stopping="private")
  out_np<-NoisyNewton(x=x,y=y,private=T,mu=1,scale=F,s0=2,beta0=rep(0,4),stepsize=0.50,maxiter=iter.j,stopping="non-private")
  
  zscores_raw[i,]<-(out_priv$beta-beta_gen)/sqrt(diag(out_priv$sandwich)/n)
  zscores_corrected[i,]<-(out_priv$beta-beta_gen)/sqrt(diag(out_priv$variances))
  zscores_np[i,]<-(out_np$beta-beta_gen)/sqrt(diag(out_np$variances))

  }

  raw_coverage_mat[j,]<-apply(zscores_raw,MARGIN = 2,coverage_calculator)
  corrected_coverage_mat[j,]<-apply(zscores_corrected,MARGIN = 2,coverage_calculator)
  np_coverage_mat[j,]<-apply(zscores_np,MARGIN = 2,coverage_calculator)
}



plot(1:7,raw_coverage_mat[,3],ylim=c(0.5,1),col=1,xaxt='n',xlab="sample size",ylab="95% CI coverage",main="",pch=16)
axis(1,at=1:7,labels=samplesizes)
abline(h=0.95,col="gray")
points(1:7,np_coverage_mat[,3],col="gray",pch=16)
points(1:7,corrected_coverage_mat[,3],col="cornflowerblue",pch=16)
legend("bottomright",legend=c("non-private stopping","private stopping, corrected","private stopping, no correction"),col=c("gray","cornflowerblue",1),pch=16,bty='n',cex=0.8)
