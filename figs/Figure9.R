source("src/extra/NGD_linear_regression_extra.R")

coverage_calculator<-function(v){
  return(mean(abs(v)<qnorm(0.975)))
}

samplesizes<-c(200,300,400,500,625,750,1000,2000,3000,4000)
niter<-c(30,33,35,41,48,51,56,60,64,65)
beta_gen=rep(1,4)
s0=2
p<-length(beta_gen)

ns<-1000

raw_coverage_mat<-matrix(NA,nrow=length(samplesizes),ncol=p)
corrected_coverage_mat<-matrix(NA,nrow=length(samplesizes),ncol=p)
np_coverage_mat<-matrix(NA,nrow=length(samplesizes),ncol=p)


set.seed(234)
for(j in 1:length(samplesizes)){
  n<-samplesizes[j]
  iter.j<-niter[j]

  zscores_raw=zscores_corrected=zscores_np=matrix(nrow=ns,ncol=p)
  for(i in 1:ns){
    x=matrix(rnorm(n*p,mean=0,sd=s0),nrow=n)
    x[,1]<-rep(1,n)
    y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
    
    out_priv<-NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=iter.j,mnorm=sqrt(2),suppress.inference=FALSE,stopping=0)
    out_np<-NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=iter.j,mnorm=sqrt(2),suppress.inference=FALSE,stopping="non-private")
  
    zscores_raw[i,]<-(out_priv$beta-beta_gen)/sqrt(diag(out_priv$sandwich)/n)
    zscores_corrected[i,]<-(out_priv$beta-beta_gen)/sqrt(diag(out_priv$variances))
    zscores_np[i,]<-(out_np$beta-beta_gen)/sqrt(diag(out_np$variances))
  }
  raw_coverage_mat[j,]<-apply(zscores_raw,MARGIN = 2,coverage_calculator)
  corrected_coverage_mat[j,]<-apply(zscores_corrected,MARGIN = 2,coverage_calculator)
  np_coverage_mat[j,]<-apply(zscores_np,MARGIN = 2,coverage_calculator)

}


plot(1:10,raw_coverage_mat[,2],xlab="sample size",ylab="95% CI coverage",xaxt="n",ylim=c(0.7,1),pch=20,col=1, main=NULL)
points(1:10,np_coverage_mat[,2],pch=20,col="gray")
points(1:10,corrected_coverage_mat[,2],pch=20,col="cornflowerblue")
axis(1, at=1:10, labels=samplesizes)
abline(h=0.95,lty=2,col="gray")
legend("bottomright",legend=c( "using true gradient", "corrected","uncorrected"),
       col=c("gray","cornflowerblue",1),pch=c(16,16,16),bty='n',cex=0.9)



