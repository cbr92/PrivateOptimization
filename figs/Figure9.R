source("NGD_linear_regression_extra.R")

coverage_calculator<-function(v){
  return(mean(abs(v)<qnorm(0.975)))
}

samplesizes<-c(200,300,400,500,625,750,1000,2000,3000,4000)
beta_gen=rep(1,4)
s0=2

ns<-1000
#niter<-c(30,33,36,42,45,52,55,57,60,63) #old
#niter<-c(30,37,41,47, , , , , , ) #1 and 2 fixed
niter<-c(30,33,35,41,48,51,56,60,60,63)
raw_coverage_mat<-matrix(NA,nrow=10,ncol=4)
corrected_coverage_mat<-matrix(NA,nrow=10,ncol=4)
np_coverage_mat<-matrix(NA,nrow=10,ncol=4)

conv_rate_priv=conv_rate_np=rep(NA,length(samplesizes))

set.seed(234)
#for(j in 1:length(samplesizes)){
for(j in 1:9){
  n<-samplesizes[j]
  iter.j<-niter[j]
  print(iter.j)
  
  zscores_raw=zscores_corrected=zscores_np=matrix(nrow=ns,ncol=4)
  conv_priv=conv_np=rep(NA,ns)
  for(i in 1:ns){
    x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
    x[,1]<-rep(1,n)
    y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
    
    #out_priv<-NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=iter.j,mnorm=sqrt(2),suppress.inference=FALSE,stopping="private")
    out_priv<-NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=iter.j,mnorm=sqrt(2),suppress.inference=FALSE,stopping=0)
    
    out_np<-NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=iter.j,mnorm=sqrt(2),suppress.inference=FALSE,stopping="non-private")
    
    conv_priv[i]<-out_priv$priv_converge
    conv_np[i]<-out_np$nonpriv_converge

    zscores_raw[i,]<-(out_priv$beta-beta_gen)/sqrt(diag(out_priv$sandwich)/n)
    zscores_corrected[i,]<-(out_priv$beta-beta_gen)/sqrt(diag(out_priv$variances))
    zscores_np[i,]<-(out_np$beta-beta_gen)/sqrt(diag(out_np$variances))
  }
  raw_coverage_mat[j,]<-apply(zscores_raw,MARGIN = 2,coverage_calculator)
  corrected_coverage_mat[j,]<-apply(zscores_corrected,MARGIN = 2,coverage_calculator)
  np_coverage_mat[j,]<-apply(zscores_np,MARGIN = 2,coverage_calculator)
  
  conv_rate_priv[j]<-mean(conv_priv)
  conv_rate_np[j]<-mean(conv_np)
}


plot(1:10,raw_coverage_mat[,2],xlab="sample size",ylab="95% CI coverage",xaxt="n",ylim=c(0.7,1),pch=20,col=1, main=NULL)
points(1:10,np_coverage_mat[,2],pch=20,col="gray")
points(1:10,corrected_coverage_mat[,2],pch=20,col="cornflowerblue")
axis(1, at=1:10, labels=samplesizes)
abline(h=0.95,lty=2,col="gray")
legend("bottomright",legend=c( "using true gradient", "corrected","uncorrected"),
       col=c("gray","cornflowerblue",1),pch=c(16,16,16),bty='n',cex=0.9)

x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
x[,1]<-rep(1,n)
y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))

out_np<-NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=30,mnorm=sqrt(2),suppress.inference=FALSE,stopping="non-private")

plot(out_np$nonpriv_gradtraj)



n=200
ns=1000

conv_yn<-rep(NA,ns)
niter<-rep(NA,ns)
betamat=varmat=zscore=alt_zscore=rawzscore_alt=rawzscore=zscore_alt1=zscore_alt2=zscore_alt3=alt_zscore_alt3=matrix(NA,ns,4)

#set.seed(170)
for(i in 1:ns){
  x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
  x[,1]<-rep(1,n)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,epsilon=8/50,delta=1/(50*n),maxiter=50)
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=50)
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=30)
  out_all=allNGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=30)
  #conv_yn[i]<-out$converge
  #niter[i]<-out$iter
  #betamat[i,]<-out$beta
  #varmat[i,1]=out$var1
  #varmat[i,2]=out$var2
  #varmat[i,3]=out$var3
  #varmat[i,4]=out$var4
  
  
  
  
  #zscore_alt1[i,1]=out$zscore1_alt1
  #zscore_alt1[i,2]=out$zscore2_alt1
  #zscore_alt1[i,3]=out$zscore3_alt1
  #zscore_alt1[i,4]=out$zscore4_alt1
  
  zscore_alt2[i,1]=out$zscore1_alt3
  zscore_alt2[i,2]=out$zscore2_alt3
  zscore_alt2[i,3]=out$zscore3_alt3
  zscore_alt2[i,4]=out$zscore4_alt3
  zscore_alt3[i,1]=out_all$zscore1_alt3
  zscore_alt3[i,2]=out_all$zscore2_alt3
  zscore_alt3[i,3]=out_all$zscore3_alt3
  zscore_alt3[i,4]=out_all$zscore4_alt3
  
  
  rawzscore[i,1]=out_all$rzscore1
  rawzscore[i,2]=out_all$rzscore2
  rawzscore[i,3]=out_all$rzscore3
  rawzscore[i,4]=out_all$rzscore4
  
  rawzscore_alt[i,1]=out$rzscore1
  rawzscore_alt[i,2]=out$rzscore2
  rawzscore_alt[i,3]=out$rzscore3
  rawzscore_alt[i,4]=out$rzscore4
  
  #out_ols<-lm(y~-1+x)
  #zscore_ols[i,]<-((summary(out_ols)$coefficients[,1])-1)/summary(out_ols)$coefficients[,2]
}
mean(abs(zscore[,2])<qnorm(0.975))


coverage_mat[1,1]<-mean(abs(zscore[,1])<qnorm(0.975))
coverage_mat[1,2]<-mean(abs(zscore[,2])<qnorm(0.975))
coverage_mat[1,3]<-mean(abs(zscore[,3])<qnorm(0.975))
coverage_mat[1,4]<-mean(abs(zscore[,4])<qnorm(0.975))

raw_coverage_mat[1,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat[1,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat[1,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat[1,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))
raw_coverage_mat_alt[1,1]<-mean(abs(rawzscore_alt[,1])<qnorm(0.975))
raw_coverage_mat_alt[1,2]<-mean(abs(rawzscore_alt[,2])<qnorm(0.975))
raw_coverage_mat_alt[1,3]<-mean(abs(rawzscore_alt[,3])<qnorm(0.975))
raw_coverage_mat_alt[1,4]<-mean(abs(rawzscore_alt[,4])<qnorm(0.975))

#alt1_coverage_mat[1,1]<-mean(abs(zscore_alt1[,1])<qnorm(0.975))
#alt1_coverage_mat[1,2]<-mean(abs(zscore_alt1[,2])<qnorm(0.975))
#alt1_coverage_mat[1,3]<-mean(abs(zscore_alt1[,3])<qnorm(0.975))
#alt1_coverage_mat[1,4]<-mean(abs(zscore_alt1[,4])<qnorm(0.975))

alt2_coverage_mat[1,1]<-mean(abs(zscore_alt2[,1])<qnorm(0.975))
alt2_coverage_mat[1,2]<-mean(abs(zscore_alt2[,2])<qnorm(0.975))
alt2_coverage_mat[1,3]<-mean(abs(zscore_alt2[,3])<qnorm(0.975))
alt2_coverage_mat[1,4]<-mean(abs(zscore_alt2[,4])<qnorm(0.975))

alt3_coverage_mat[1,1]<-mean(abs(zscore_alt3[,1])<qnorm(0.975))
alt3_coverage_mat[1,2]<-mean(abs(zscore_alt3[,2])<qnorm(0.975))
alt3_coverage_mat[1,3]<-mean(abs(zscore_alt3[,3])<qnorm(0.975))
alt3_coverage_mat[1,4]<-mean(abs(zscore_alt3[,4])<qnorm(0.975))


n=300
ns=1000

conv_yn<-rep(NA,ns)
niter<-rep(NA,ns)
betamat=varmat=zscore=rawzscore=zscore_alt1=zscore_alt2=zscore_alt3=matrix(NA,ns,4)

#set.seed(278)
for(i in 1:ns){
  x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
  x[,1]<-rep(1,n)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,epsilon=8/50,delta=1/(50*n),maxiter=50)
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=50)
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=33)
  conv_yn[i]<-out$converge
  niter[i]<-out$iter
  betamat[i,]<-out$beta
  #varmat[i,1]=out$var1
  #varmat[i,2]=out$var2
  #varmat[i,3]=out$var3
  #varmat[i,4]=out$var4
  
  zscore[i,1]=out$zscore1
  zscore[i,2]=out$zscore2
  zscore[i,3]=out$zscore3
  zscore[i,4]=out$zscore4
  
  zscore_alt1[i,1]=out$zscore1_alt1
  zscore_alt1[i,2]=out$zscore2_alt1
  zscore_alt1[i,3]=out$zscore3_alt1
  zscore_alt1[i,4]=out$zscore4_alt1
  
  zscore_alt2[i,1]=out$zscore1_alt2
  zscore_alt2[i,2]=out$zscore2_alt2
  zscore_alt2[i,3]=out$zscore3_alt2
  zscore_alt2[i,4]=out$zscore4_alt2
  zscore_alt3[i,1]=out$zscore1_alt3
  zscore_alt3[i,2]=out$zscore2_alt3
  zscore_alt3[i,3]=out$zscore3_alt3
  zscore_alt3[i,4]=out$zscore4_alt3
  
  rawzscore[i,1]=out$rzscore1
  rawzscore[i,2]=out$rzscore2
  rawzscore[i,3]=out$rzscore3
  rawzscore[i,4]=out$rzscore4
  #out_ols<-lm(y~-1+x)
  #zscore_ols[i,]<-((summary(out_ols)$coefficients[,1])-1)/summary(out_ols)$coefficients[,2]
}
mean(abs(zscore[,2])<qnorm(0.975))


coverage_mat[2,1]<-mean(abs(zscore[,1])<qnorm(0.975))
coverage_mat[2,2]<-mean(abs(zscore[,2])<qnorm(0.975))
coverage_mat[2,3]<-mean(abs(zscore[,3])<qnorm(0.975))
coverage_mat[2,4]<-mean(abs(zscore[,4])<qnorm(0.975))
raw_coverage_mat[2,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat[2,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat[2,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat[2,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))
raw_coverage_mat_alt[2,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat_alt[2,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat_alt[2,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat_alt[2,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))

alt1_coverage_mat[2,1]<-mean(abs(zscore_alt1[,1])<qnorm(0.975))
alt1_coverage_mat[2,2]<-mean(abs(zscore_alt1[,2])<qnorm(0.975))
alt1_coverage_mat[2,3]<-mean(abs(zscore_alt1[,3])<qnorm(0.975))
alt1_coverage_mat[2,4]<-mean(abs(zscore_alt1[,4])<qnorm(0.975))

alt2_coverage_mat[2,1]<-mean(abs(zscore_alt2[,1])<qnorm(0.975))
alt2_coverage_mat[2,2]<-mean(abs(zscore_alt2[,2])<qnorm(0.975))
alt2_coverage_mat[2,3]<-mean(abs(zscore_alt2[,3])<qnorm(0.975))
alt2_coverage_mat[2,4]<-mean(abs(zscore_alt2[,4])<qnorm(0.975))
alt3_coverage_mat[2,1]<-mean(abs(zscore_alt3[,1])<qnorm(0.975))
alt3_coverage_mat[2,2]<-mean(abs(zscore_alt3[,2])<qnorm(0.975))
alt3_coverage_mat[2,3]<-mean(abs(zscore_alt3[,3])<qnorm(0.975))
alt3_coverage_mat[2,4]<-mean(abs(zscore_alt3[,4])<qnorm(0.975))


n=400
ns=1000

conv_yn<-rep(NA,ns)
niter<-rep(NA,ns)
betamat=varmat=zscore=rawzscore=zscore_ols=matrix(NA,ns,4)

#set.seed(914)
for(i in 1:ns){
  x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
  x[,1]<-rep(1,n)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,epsilon=8/50,delta=1/(50*n),maxiter=50)
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=50)
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=36)
  conv_yn[i]<-out$converge
  niter[i]<-out$iter
  betamat[i,]<-out$beta
  #varmat[i,1]=out$var1
  #varmat[i,2]=out$var2
  #varmat[i,3]=out$var3
  #varmat[i,4]=out$var4
  
  zscore[i,1]=out$zscore1
  zscore[i,2]=out$zscore2
  zscore[i,3]=out$zscore3
  zscore[i,4]=out$zscore4
  
  zscore_alt1[i,1]=out$zscore1_alt1
  zscore_alt1[i,2]=out$zscore2_alt1
  zscore_alt1[i,3]=out$zscore3_alt1
  zscore_alt1[i,4]=out$zscore4_alt1
  
  zscore_alt2[i,1]=out$zscore1_alt2
  zscore_alt2[i,2]=out$zscore2_alt2
  zscore_alt2[i,3]=out$zscore3_alt2
  zscore_alt2[i,4]=out$zscore4_alt2
  zscore_alt3[i,1]=out$zscore1_alt3
  zscore_alt3[i,2]=out$zscore2_alt3
  zscore_alt3[i,3]=out$zscore3_alt3
  zscore_alt3[i,4]=out$zscore4_alt3
  
  rawzscore[i,1]=out$rzscore1
  rawzscore[i,2]=out$rzscore2
  rawzscore[i,3]=out$rzscore3
  rawzscore[i,4]=out$rzscore4
  #out_ols<-lm(y~-1+x)
  #zscore_ols[i,]<-((summary(out_ols)$coefficients[,1])-1)/summary(out_ols)$coefficients[,2]
}
mean(abs(zscore[,2])<qnorm(0.975))


coverage_mat[3,1]<-mean(abs(zscore[,1])<qnorm(0.975))
coverage_mat[3,2]<-mean(abs(zscore[,2])<qnorm(0.975))
coverage_mat[3,3]<-mean(abs(zscore[,3])<qnorm(0.975))
coverage_mat[3,4]<-mean(abs(zscore[,4])<qnorm(0.975))
raw_coverage_mat[3,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat[3,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat[3,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat[3,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))

alt1_coverage_mat[3,1]<-mean(abs(zscore_alt1[,1])<qnorm(0.975))
alt1_coverage_mat[3,2]<-mean(abs(zscore_alt1[,2])<qnorm(0.975))
alt1_coverage_mat[3,3]<-mean(abs(zscore_alt1[,3])<qnorm(0.975))
alt1_coverage_mat[3,4]<-mean(abs(zscore_alt1[,4])<qnorm(0.975))

alt2_coverage_mat[3,1]<-mean(abs(zscore_alt2[,1])<qnorm(0.975))
alt2_coverage_mat[3,2]<-mean(abs(zscore_alt2[,2])<qnorm(0.975))
alt2_coverage_mat[3,3]<-mean(abs(zscore_alt2[,3])<qnorm(0.975))
alt2_coverage_mat[3,4]<-mean(abs(zscore_alt2[,4])<qnorm(0.975))
alt3_coverage_mat[3,1]<-mean(abs(zscore_alt3[,1])<qnorm(0.975))
alt3_coverage_mat[3,2]<-mean(abs(zscore_alt3[,2])<qnorm(0.975))
alt3_coverage_mat[3,3]<-mean(abs(zscore_alt3[,3])<qnorm(0.975))
alt3_coverage_mat[3,4]<-mean(abs(zscore_alt3[,4])<qnorm(0.975))


n=500
ns=1000

conv_yn<-rep(NA,ns)
niter<-rep(NA,ns)
betamat=varmat=zscore=rawzscore=zscore_ols=matrix(NA,ns,4)

#set.seed(1002)
for(i in 1:ns){
  x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
  x[,1]<-rep(1,n)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,epsilon=8/50,delta=1/(50*n),maxiter=50)
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=42)
  conv_yn[i]<-out$converge
  niter[i]<-out$iter
  betamat[i,]<-out$beta
  #varmat[i,1]=out$var1
  #varmat[i,2]=out$var2
  #varmat[i,3]=out$var3
  #varmat[i,4]=out$var4
  
  zscore[i,1]=out$zscore1
  zscore[i,2]=out$zscore2
  zscore[i,3]=out$zscore3
  zscore[i,4]=out$zscore4
  
  zscore_alt1[i,1]=out$zscore1_alt1
  zscore_alt1[i,2]=out$zscore2_alt1
  zscore_alt1[i,3]=out$zscore3_alt1
  zscore_alt1[i,4]=out$zscore4_alt1
  
  zscore_alt2[i,1]=out$zscore1_alt2
  zscore_alt2[i,2]=out$zscore2_alt2
  zscore_alt2[i,3]=out$zscore3_alt2
  zscore_alt2[i,4]=out$zscore4_alt2
  zscore_alt3[i,1]=out$zscore1_alt3
  zscore_alt3[i,2]=out$zscore2_alt3
  zscore_alt3[i,3]=out$zscore3_alt3
  zscore_alt3[i,4]=out$zscore4_alt3
  
  rawzscore[i,1]=out$rzscore1
  rawzscore[i,2]=out$rzscore2
  rawzscore[i,3]=out$rzscore3
  rawzscore[i,4]=out$rzscore4
  #out_ols<-lm(y~-1+x)
  #zscore_ols[i,]<-((summary(out_ols)$coefficients[,1])-1)/summary(out_ols)$coefficients[,2]
}
mean(abs(zscore[,2])<qnorm(0.975))


coverage_mat[4,1]<-mean(abs(zscore[,1])<qnorm(0.975))
coverage_mat[4,2]<-mean(abs(zscore[,2])<qnorm(0.975))
coverage_mat[4,3]<-mean(abs(zscore[,3])<qnorm(0.975))
coverage_mat[4,4]<-mean(abs(zscore[,4])<qnorm(0.975))
raw_coverage_mat[4,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat[4,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat[4,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat[4,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))

alt1_coverage_mat[4,1]<-mean(abs(zscore_alt1[,1])<qnorm(0.975))
alt1_coverage_mat[4,2]<-mean(abs(zscore_alt1[,2])<qnorm(0.975))
alt1_coverage_mat[4,3]<-mean(abs(zscore_alt1[,3])<qnorm(0.975))
alt1_coverage_mat[4,4]<-mean(abs(zscore_alt1[,4])<qnorm(0.975))

alt2_coverage_mat[4,1]<-mean(abs(zscore_alt2[,1])<qnorm(0.975))
alt2_coverage_mat[4,2]<-mean(abs(zscore_alt2[,2])<qnorm(0.975))
alt2_coverage_mat[4,3]<-mean(abs(zscore_alt2[,3])<qnorm(0.975))
alt2_coverage_mat[4,4]<-mean(abs(zscore_alt2[,4])<qnorm(0.975))
alt3_coverage_mat[4,1]<-mean(abs(zscore_alt3[,1])<qnorm(0.975))
alt3_coverage_mat[4,2]<-mean(abs(zscore_alt3[,2])<qnorm(0.975))
alt3_coverage_mat[4,3]<-mean(abs(zscore_alt3[,3])<qnorm(0.975))
alt3_coverage_mat[4,4]<-mean(abs(zscore_alt3[,4])<qnorm(0.975))




n=625
ns=1000

conv_yn<-rep(NA,ns)
niter<-rep(NA,ns)
betamat=varmat=zscore=rawzscore=zscore_ols=matrix(NA,ns,4)

#set.seed(520)
for(i in 1:ns){
  x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
  x[,1]<-rep(1,n)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,epsilon=8/50,delta=1/(50*n),maxiter=50)
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=50)
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=45)
  conv_yn[i]<-out$converge
  niter[i]<-out$iter
  betamat[i,]<-out$beta
  #varmat[i,1]=out$var1
  #varmat[i,2]=out$var2
  #varmat[i,3]=out$var3
  #varmat[i,4]=out$var4
  
  zscore[i,1]=out$zscore1
  zscore[i,2]=out$zscore2
  zscore[i,3]=out$zscore3
  zscore[i,4]=out$zscore4
  
  zscore_alt1[i,1]=out$zscore1_alt1
  zscore_alt1[i,2]=out$zscore2_alt1
  zscore_alt1[i,3]=out$zscore3_alt1
  zscore_alt1[i,4]=out$zscore4_alt1
  
  zscore_alt2[i,1]=out$zscore1_alt2
  zscore_alt2[i,2]=out$zscore2_alt2
  zscore_alt2[i,3]=out$zscore3_alt2
  zscore_alt2[i,4]=out$zscore4_alt2
  zscore_alt3[i,1]=out$zscore1_alt3
  zscore_alt3[i,2]=out$zscore2_alt3
  zscore_alt3[i,3]=out$zscore3_alt3
  zscore_alt3[i,4]=out$zscore4_alt3
  
  rawzscore[i,1]=out$rzscore1
  rawzscore[i,2]=out$rzscore2
  rawzscore[i,3]=out$rzscore3
  rawzscore[i,4]=out$rzscore4
  #out_ols<-lm(y~-1+x)
  #zscore_ols[i,]<-((summary(out_ols)$coefficients[,1])-1)/summary(out_ols)$coefficients[,2]
}
mean(abs(zscore[,2])<qnorm(0.975))


coverage_mat[5,1]<-mean(abs(zscore[,1])<qnorm(0.975))
coverage_mat[5,2]<-mean(abs(zscore[,2])<qnorm(0.975))
coverage_mat[5,3]<-mean(abs(zscore[,3])<qnorm(0.975))
coverage_mat[5,4]<-mean(abs(zscore[,4])<qnorm(0.975))
raw_coverage_mat[5,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat[5,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat[5,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat[5,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))

alt1_coverage_mat[5,1]<-mean(abs(zscore_alt1[,1])<qnorm(0.975))
alt1_coverage_mat[5,2]<-mean(abs(zscore_alt1[,2])<qnorm(0.975))
alt1_coverage_mat[5,3]<-mean(abs(zscore_alt1[,3])<qnorm(0.975))
alt1_coverage_mat[5,4]<-mean(abs(zscore_alt1[,4])<qnorm(0.975))

alt2_coverage_mat[5,1]<-mean(abs(zscore_alt2[,1])<qnorm(0.975))
alt2_coverage_mat[5,2]<-mean(abs(zscore_alt2[,2])<qnorm(0.975))
alt2_coverage_mat[5,3]<-mean(abs(zscore_alt2[,3])<qnorm(0.975))
alt2_coverage_mat[5,4]<-mean(abs(zscore_alt2[,4])<qnorm(0.975))
alt3_coverage_mat[5,1]<-mean(abs(zscore_alt3[,1])<qnorm(0.975))
alt3_coverage_mat[5,2]<-mean(abs(zscore_alt3[,2])<qnorm(0.975))
alt3_coverage_mat[5,3]<-mean(abs(zscore_alt3[,3])<qnorm(0.975))
alt3_coverage_mat[5,4]<-mean(abs(zscore_alt3[,4])<qnorm(0.975))



n=750
ns=1000

conv_yn<-rep(NA,ns)
niter<-rep(NA,ns)
betamat=varmat=zscore=rawzscore=zscore_ols=matrix(NA,ns,4)

#set.seed(1225)
for(i in 1:ns){
  x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
  x[,1]<-rep(1,n)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,epsilon=8/50,delta=1/(50*n),maxiter=50)
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=50)
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=52)
  conv_yn[i]<-out$converge
  niter[i]<-out$iter
  betamat[i,]<-out$beta
  #varmat[i,1]=out$var1
  #varmat[i,2]=out$var2
  #varmat[i,3]=out$var3
  #varmat[i,4]=out$var4
  
  zscore[i,1]=out$zscore1
  zscore[i,2]=out$zscore2
  zscore[i,3]=out$zscore3
  zscore[i,4]=out$zscore4
  
  zscore_alt1[i,1]=out$zscore1_alt1
  zscore_alt1[i,2]=out$zscore2_alt1
  zscore_alt1[i,3]=out$zscore3_alt1
  zscore_alt1[i,4]=out$zscore4_alt1
  
  zscore_alt2[i,1]=out$zscore1_alt2
  zscore_alt2[i,2]=out$zscore2_alt2
  zscore_alt2[i,3]=out$zscore3_alt2
  zscore_alt2[i,4]=out$zscore4_alt2
  zscore_alt3[i,1]=out$zscore1_alt3
  zscore_alt3[i,2]=out$zscore2_alt3
  zscore_alt3[i,3]=out$zscore3_alt3
  zscore_alt3[i,4]=out$zscore4_alt3
  
  rawzscore[i,1]=out$rzscore1
  rawzscore[i,2]=out$rzscore2
  rawzscore[i,3]=out$rzscore3
  rawzscore[i,4]=out$rzscore4
  #out_ols<-lm(y~-1+x)
  #zscore_ols[i,]<-((summary(out_ols)$coefficients[,1])-1)/summary(out_ols)$coefficients[,2]
}
mean(abs(zscore[,2])<qnorm(0.975))


coverage_mat[6,1]<-mean(abs(zscore[,1])<qnorm(0.975))
coverage_mat[6,2]<-mean(abs(zscore[,2])<qnorm(0.975))
coverage_mat[6,3]<-mean(abs(zscore[,3])<qnorm(0.975))
coverage_mat[6,4]<-mean(abs(zscore[,4])<qnorm(0.975))
raw_coverage_mat[6,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat[6,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat[6,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat[6,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))

alt1_coverage_mat[6,1]<-mean(abs(zscore_alt1[,1])<qnorm(0.975))
alt1_coverage_mat[6,2]<-mean(abs(zscore_alt1[,2])<qnorm(0.975))
alt1_coverage_mat[6,3]<-mean(abs(zscore_alt1[,3])<qnorm(0.975))
alt1_coverage_mat[6,4]<-mean(abs(zscore_alt1[,4])<qnorm(0.975))

alt2_coverage_mat[6,1]<-mean(abs(zscore_alt2[,1])<qnorm(0.975))
alt2_coverage_mat[6,2]<-mean(abs(zscore_alt2[,2])<qnorm(0.975))
alt2_coverage_mat[6,3]<-mean(abs(zscore_alt2[,3])<qnorm(0.975))
alt2_coverage_mat[6,4]<-mean(abs(zscore_alt2[,4])<qnorm(0.975))
alt3_coverage_mat[6,1]<-mean(abs(zscore_alt3[,1])<qnorm(0.975))
alt3_coverage_mat[6,2]<-mean(abs(zscore_alt3[,2])<qnorm(0.975))
alt3_coverage_mat[6,3]<-mean(abs(zscore_alt3[,3])<qnorm(0.975))
alt3_coverage_mat[6,4]<-mean(abs(zscore_alt3[,4])<qnorm(0.975))



n=1000
ns=1000

conv_yn<-rep(NA,ns)
niter<-rep(NA,ns)
betamat=varmat=zscore=rawzscore=zscore_ols=matrix(NA,ns,4)

#set.seed(517)
for(i in 1:ns){
  x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
  x[,1]<-rep(1,n)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,epsilon=8/50,delta=1/(50*n),maxiter=50)
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=50)
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=55)
  conv_yn[i]<-out$converge
  niter[i]<-out$iter
  betamat[i,]<-out$beta
  #varmat[i,1]=out$var1
  #varmat[i,2]=out$var2
  #varmat[i,3]=out$var3
  #varmat[i,4]=out$var4
  
  zscore[i,1]=out$zscore1
  zscore[i,2]=out$zscore2
  zscore[i,3]=out$zscore3
  zscore[i,4]=out$zscore4
  
  zscore_alt1[i,1]=out$zscore1_alt1
  zscore_alt1[i,2]=out$zscore2_alt1
  zscore_alt1[i,3]=out$zscore3_alt1
  zscore_alt1[i,4]=out$zscore4_alt1
  
  zscore_alt2[i,1]=out$zscore1_alt2
  zscore_alt2[i,2]=out$zscore2_alt2
  zscore_alt2[i,3]=out$zscore3_alt2
  zscore_alt2[i,4]=out$zscore4_alt2
  zscore_alt3[i,1]=out$zscore1_alt3
  zscore_alt3[i,2]=out$zscore2_alt3
  zscore_alt3[i,3]=out$zscore3_alt3
  zscore_alt3[i,4]=out$zscore4_alt3
  
  rawzscore[i,1]=out$rzscore1
  rawzscore[i,2]=out$rzscore2
  rawzscore[i,3]=out$rzscore3
  rawzscore[i,4]=out$rzscore4
  #out_ols<-lm(y~-1+x)
  #zscore_ols[i,]<-((summary(out_ols)$coefficients[,1])-1)/summary(out_ols)$coefficients[,2]
}
mean(abs(zscore[,2])<qnorm(0.975))


coverage_mat[7,1]<-mean(abs(zscore[,1])<qnorm(0.975))
coverage_mat[7,2]<-mean(abs(zscore[,2])<qnorm(0.975))
coverage_mat[7,3]<-mean(abs(zscore[,3])<qnorm(0.975))
coverage_mat[7,4]<-mean(abs(zscore[,4])<qnorm(0.975))
raw_coverage_mat[7,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat[7,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat[7,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat[7,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))

alt1_coverage_mat[7,1]<-mean(abs(zscore_alt1[,1])<qnorm(0.975))
alt1_coverage_mat[7,2]<-mean(abs(zscore_alt1[,2])<qnorm(0.975))
alt1_coverage_mat[7,3]<-mean(abs(zscore_alt1[,3])<qnorm(0.975))
alt1_coverage_mat[7,4]<-mean(abs(zscore_alt1[,4])<qnorm(0.975))

alt2_coverage_mat[7,1]<-mean(abs(zscore_alt2[,1])<qnorm(0.975))
alt2_coverage_mat[7,2]<-mean(abs(zscore_alt2[,2])<qnorm(0.975))
alt2_coverage_mat[7,3]<-mean(abs(zscore_alt2[,3])<qnorm(0.975))
alt2_coverage_mat[7,4]<-mean(abs(zscore_alt2[,4])<qnorm(0.975))
alt3_coverage_mat[7,1]<-mean(abs(zscore_alt3[,1])<qnorm(0.975))
alt3_coverage_mat[7,2]<-mean(abs(zscore_alt3[,2])<qnorm(0.975))
alt3_coverage_mat[7,3]<-mean(abs(zscore_alt3[,3])<qnorm(0.975))
alt3_coverage_mat[7,4]<-mean(abs(zscore_alt3[,4])<qnorm(0.975))




n=2000
ns=1000

conv_yn<-rep(NA,ns)
niter<-rep(NA,ns)
betamat=varmat=zscore=rawzscore=zscore_ols=matrix(NA,ns,4)

#set.seed(641)
for(i in 1:ns){
  x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
  x[,1]<-rep(1,n)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,epsilon=8/50,delta=1/(50*n),maxiter=50)
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=57)
  conv_yn[i]<-out$converge
  niter[i]<-out$iter
  betamat[i,]<-out$beta
  #varmat[i,1]=out$var1
  #varmat[i,2]=out$var2
  #varmat[i,3]=out$var3
  #varmat[i,4]=out$var4
  
  zscore[i,1]=out$zscore1
  zscore[i,2]=out$zscore2
  zscore[i,3]=out$zscore3
  zscore[i,4]=out$zscore4
  
  zscore_alt1[i,1]=out$zscore1_alt1
  zscore_alt1[i,2]=out$zscore2_alt1
  zscore_alt1[i,3]=out$zscore3_alt1
  zscore_alt1[i,4]=out$zscore4_alt1
  
  zscore_alt2[i,1]=out$zscore1_alt2
  zscore_alt2[i,2]=out$zscore2_alt2
  zscore_alt2[i,3]=out$zscore3_alt2
  zscore_alt2[i,4]=out$zscore4_alt2
  zscore_alt3[i,1]=out$zscore1_alt3
  zscore_alt3[i,2]=out$zscore2_alt3
  zscore_alt3[i,3]=out$zscore3_alt3
  zscore_alt3[i,4]=out$zscore4_alt3
  
  rawzscore[i,1]=out$rzscore1
  rawzscore[i,2]=out$rzscore2
  rawzscore[i,3]=out$rzscore3
  rawzscore[i,4]=out$rzscore4
  #out_ols<-lm(y~-1+x)
  #zscore_ols[i,]<-((summary(out_ols)$coefficients[,1])-1)/summary(out_ols)$coefficients[,2]
}
mean(abs(zscore[,2])<qnorm(0.975))


coverage_mat[8,1]<-mean(abs(zscore[,1])<qnorm(0.975))
coverage_mat[8,2]<-mean(abs(zscore[,2])<qnorm(0.975))
coverage_mat[8,3]<-mean(abs(zscore[,3])<qnorm(0.975))
coverage_mat[8,4]<-mean(abs(zscore[,4])<qnorm(0.975))
raw_coverage_mat[8,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat[8,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat[8,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat[8,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))

alt1_coverage_mat[8,1]<-mean(abs(zscore_alt1[,1])<qnorm(0.975))
alt1_coverage_mat[8,2]<-mean(abs(zscore_alt1[,2])<qnorm(0.975))
alt1_coverage_mat[8,3]<-mean(abs(zscore_alt1[,3])<qnorm(0.975))
alt1_coverage_mat[8,4]<-mean(abs(zscore_alt1[,4])<qnorm(0.975))

alt2_coverage_mat[8,1]<-mean(abs(zscore_alt2[,1])<qnorm(0.975))
alt2_coverage_mat[8,2]<-mean(abs(zscore_alt2[,2])<qnorm(0.975))
alt2_coverage_mat[8,3]<-mean(abs(zscore_alt2[,3])<qnorm(0.975))
alt2_coverage_mat[8,4]<-mean(abs(zscore_alt2[,4])<qnorm(0.975))
alt3_coverage_mat[8,1]<-mean(abs(zscore_alt3[,1])<qnorm(0.975))
alt3_coverage_mat[8,2]<-mean(abs(zscore_alt3[,2])<qnorm(0.975))
alt3_coverage_mat[8,3]<-mean(abs(zscore_alt3[,3])<qnorm(0.975))
alt3_coverage_mat[8,4]<-mean(abs(zscore_alt3[,4])<qnorm(0.975))




n=3000
ns=1000

conv_yn<-rep(NA,ns)
niter<-rep(NA,ns)
betamat=varmat=zscore=rawzscore=zscore_ols=matrix(NA,ns,4)

#set.seed(833)
for(i in 1:ns){
  x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
  x[,1]<-rep(1,n)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,epsilon=8/50,delta=1/(50*n),maxiter=50)
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=75)
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=60)
  conv_yn[i]<-out$converge
  niter[i]<-out$iter
  betamat[i,]<-out$beta
  #varmat[i,1]=out$var1
  #varmat[i,2]=out$var2
  #varmat[i,3]=out$var3
  #varmat[i,4]=out$var4
  
  zscore[i,1]=out$zscore1
  zscore[i,2]=out$zscore2
  zscore[i,3]=out$zscore3
  zscore[i,4]=out$zscore4
  
  zscore_alt1[i,1]=out$zscore1_alt1
  zscore_alt1[i,2]=out$zscore2_alt1
  zscore_alt1[i,3]=out$zscore3_alt1
  zscore_alt1[i,4]=out$zscore4_alt1
  
  zscore_alt2[i,1]=out$zscore1_alt2
  zscore_alt2[i,2]=out$zscore2_alt2
  zscore_alt2[i,3]=out$zscore3_alt2
  zscore_alt2[i,4]=out$zscore4_alt2
  zscore_alt3[i,1]=out$zscore1_alt3
  zscore_alt3[i,2]=out$zscore2_alt3
  zscore_alt3[i,3]=out$zscore3_alt3
  zscore_alt3[i,4]=out$zscore4_alt3
  
  rawzscore[i,1]=out$rzscore1
  rawzscore[i,2]=out$rzscore2
  rawzscore[i,3]=out$rzscore3
  rawzscore[i,4]=out$rzscore4
  #out_ols<-lm(y~-1+x)
  #zscore_ols[i,]<-((summary(out_ols)$coefficients[,1])-1)/summary(out_ols)$coefficients[,2]
}
mean(abs(zscore[,2])<qnorm(0.975))


coverage_mat[9,1]<-mean(abs(zscore[,1])<qnorm(0.975))
coverage_mat[9,2]<-mean(abs(zscore[,2])<qnorm(0.975))
coverage_mat[9,3]<-mean(abs(zscore[,3])<qnorm(0.975))
coverage_mat[9,4]<-mean(abs(zscore[,4])<qnorm(0.975))
raw_coverage_mat[9,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat[9,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat[9,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat[9,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))

alt1_coverage_mat[9,1]<-mean(abs(zscore_alt1[,1])<qnorm(0.975))
alt1_coverage_mat[9,2]<-mean(abs(zscore_alt1[,2])<qnorm(0.975))
alt1_coverage_mat[9,3]<-mean(abs(zscore_alt1[,3])<qnorm(0.975))
alt1_coverage_mat[9,4]<-mean(abs(zscore_alt1[,4])<qnorm(0.975))

alt2_coverage_mat[9,1]<-mean(abs(zscore_alt2[,1])<qnorm(0.975))
alt2_coverage_mat[9,2]<-mean(abs(zscore_alt2[,2])<qnorm(0.975))
alt2_coverage_mat[9,3]<-mean(abs(zscore_alt2[,3])<qnorm(0.975))
alt2_coverage_mat[9,4]<-mean(abs(zscore_alt2[,4])<qnorm(0.975))
alt3_coverage_mat[9,1]<-mean(abs(zscore_alt3[,1])<qnorm(0.975))
alt3_coverage_mat[9,2]<-mean(abs(zscore_alt3[,2])<qnorm(0.975))
alt3_coverage_mat[9,3]<-mean(abs(zscore_alt3[,3])<qnorm(0.975))
alt3_coverage_mat[9,4]<-mean(abs(zscore_alt3[,4])<qnorm(0.975))



n=4000
ns=1000

conv_yn<-rep(NA,ns)
niter<-rep(NA,ns)
betamat=varmat=zscore=rawzscore=zscore_ols=zscore_alt1=zscore_alt2=zscore_alt3=matrix(NA,ns,4)

#set.seed(5309)
for(i in 1:ns){
  x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
  x[,1]<-rep(1,n)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,epsilon=8/50,delta=1/(50*n),maxiter=50)
  #out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=75)
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=63)
  conv_yn[i]<-out$converge
  niter[i]<-out$iter
  betamat[i,]<-out$beta
  #varmat[i,1]=out$var1
  #varmat[i,2]=out$var2
  #varmat[i,3]=out$var3
  #varmat[i,4]=out$var4
  
  zscore[i,1]=out$zscore1
  zscore[i,2]=out$zscore2
  zscore[i,3]=out$zscore3
  zscore[i,4]=out$zscore4
  
  zscore_alt1[i,1]=out$zscore1_alt1
  zscore_alt1[i,2]=out$zscore2_alt1
  zscore_alt1[i,3]=out$zscore3_alt1
  zscore_alt1[i,4]=out$zscore4_alt1
  
  zscore_alt2[i,1]=out$zscore1_alt2
  zscore_alt2[i,2]=out$zscore2_alt2
  zscore_alt2[i,3]=out$zscore3_alt2
  zscore_alt2[i,4]=out$zscore4_alt2
  zscore_alt3[i,1]=out$zscore1_alt3
  zscore_alt3[i,2]=out$zscore2_alt3
  zscore_alt3[i,3]=out$zscore3_alt3
  zscore_alt3[i,4]=out$zscore4_alt3
  
  rawzscore[i,1]=out$rzscore1
  rawzscore[i,2]=out$rzscore2
  rawzscore[i,3]=out$rzscore3
  rawzscore[i,4]=out$rzscore4
  #out_ols<-lm(y~-1+x)
  #zscore_ols[i,]<-((summary(out_ols)$coefficients[,1])-1)/summary(out_ols)$coefficients[,2]
}
mean(abs(zscore[,2])<qnorm(0.975))


coverage_mat[10,1]<-mean(abs(zscore[,1])<qnorm(0.975))
coverage_mat[10,2]<-mean(abs(zscore[,2])<qnorm(0.975))
coverage_mat[10,3]<-mean(abs(zscore[,3])<qnorm(0.975))
coverage_mat[10,4]<-mean(abs(zscore[,4])<qnorm(0.975))
raw_coverage_mat[10,1]<-mean(abs(rawzscore[,1])<qnorm(0.975))
raw_coverage_mat[10,2]<-mean(abs(rawzscore[,2])<qnorm(0.975))
raw_coverage_mat[10,3]<-mean(abs(rawzscore[,3])<qnorm(0.975))
raw_coverage_mat[10,4]<-mean(abs(rawzscore[,4])<qnorm(0.975))

alt1_coverage_mat[10,1]<-mean(abs(zscore_alt1[,1])<qnorm(0.975))
alt1_coverage_mat[10,2]<-mean(abs(zscore_alt1[,2])<qnorm(0.975))
alt1_coverage_mat[10,3]<-mean(abs(zscore_alt1[,3])<qnorm(0.975))
alt1_coverage_mat[10,4]<-mean(abs(zscore_alt1[,4])<qnorm(0.975))

alt2_coverage_mat[10,1]<-mean(abs(zscore_alt2[,1])<qnorm(0.975))
alt2_coverage_mat[10,2]<-mean(abs(zscore_alt2[,2])<qnorm(0.975))
alt2_coverage_mat[10,3]<-mean(abs(zscore_alt2[,3])<qnorm(0.975))
alt2_coverage_mat[10,4]<-mean(abs(zscore_alt2[,4])<qnorm(0.975))
alt3_coverage_mat[10,1]<-mean(abs(zscore_alt3[,1])<qnorm(0.975))
alt3_coverage_mat[10,2]<-mean(abs(zscore_alt3[,2])<qnorm(0.975))
alt3_coverage_mat[10,3]<-mean(abs(zscore_alt3[,3])<qnorm(0.975))
alt3_coverage_mat[10,4]<-mean(abs(zscore_alt3[,4])<qnorm(0.975))



samplesizes2<-c(200,300,400,500,625,750,1000,2000,3000,4000)



save(samplesizes2,raw_coverage_mat,alt2_coverage_mat,alt3_coverage_mat,beta_gen,s0,ns,file="/Users/caseybradshaw/Documents/privacy/march_coverage_data.Rdata")

plot(1:10,raw_coverage_mat[,4],xlab="sample size",ylab="95% CI coverage",xaxt="n",ylim=c(0.7,1),pch=20,col=1, main=NULL)
points(1:10,alt2_coverage_mat[,4],pch=20,col="gray")
points(1:10,alt3_coverage_mat[,4],pch=20,col="cornflowerblue")
axis(1, at=1:10, labels=samplesizes2)
abline(h=0.95,lty=2,col="gray")
legend("bottomright",legend=c( "using true gradient", "corrected","uncorrected"),
       col=c("gray","cornflowerblue",1),pch=c(16,16,16),bty='n',cex=0.9)