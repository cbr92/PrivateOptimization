#############################################
########## CONTENTS #########################
#############################################
## 1. generate figure [2]
## 2. generate figure [3]
#############################################
#############################################

source("src/extra/clipped_gradient_descent.R")


#############################################
#########  1 (linear regression)  ###########
#############################################


source("src/NGD_linear_regression.R")

beta_gen=c(1.5,1,-1,0.5)
ns=400
N=c(250,500,1000,2000,4000)

betas_clip=betas_ngd=vector("list",length=length(N))

set.seed(1020)
for(j in 1:5)
{
  betas_clip[[j]]<-matrix(nrow=ns,ncol=4)
  betas_ngd[[j]]<-matrix(nrow=ns,ncol=4)
  
  n=N[j]
  iter.j<-(50+25*(n>500)+25*(n>2000))
  for(i in 1:ns){
    x=matrix(rnorm(n*4,mean=0,sd=1),nrow=n)
    x[,1]<-rep(1,n)
    y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=2))
    out_ngd<-NGD.Huber(x=x,y=y,private=T,scale=F,s0=2,maxiter=iter.j,mu=(0.5*sqrt((iter.j+2)/iter.j)),
                           stepsize=(1/(2*1.345))-0.001,stopping=0,mnorm=sqrt(2),suppress.inference=TRUE)
    #privacy budget in the above is calibrated to make the estimation process 0.5-GDP  (without allocating privacy
    #budget to inference, to be comparable to the clipping version)
    out_clip<-clipped_linear(x=x,y=y,private=T,scale=F,s0=2,maxiter=iter.j,mu=0.5,c_grad=1,
                          stepsize=(1/(2*1.345))-0.001,stopping=0)
    
    betas_ngd[[j]][i,]<-out_ngd$beta
    betas_clip[[j]][i,]<-out_clip$beta
  }
}

#############################################
# to generate Figure [2]
#############################################
par(mfrow=c(1,2))
boxplot(betas_ngd[[1]][,2],betas_ngd[[2]][,2],betas_ngd[[3]][,2],betas_ngd[[4]][,2],betas_ngd[[5]][,2],
        names=N, xlab="sample size",ylab=expression(paste("Estimate of ", beta[2]," parameter")), 
        main="Linear Regression, Mallows Weights",ylim=c(0.2,1.7))
abline(h=beta_gen[2],col=2)
boxplot(betas_clip[[1]][,2],betas_clip[[2]][,2],betas_clip[[3]][,2],betas_clip[[4]][,2], betas_clip[[5]][,2],
        names=N, xlab="sample size",ylab=expression(paste("Estimate of ", beta[2]," parameter")), 
        main="Linear Regression, Gradient Clipping",ylim=c(0.2,1.7))
abline(h=beta_gen[2],col=2)



#############################################
######### 2 (logistic regression) ###########
#############################################


source("src/private_logistic_regression.R")


beta_gen=c(1.5,1,-1,0.5)
ns=400
N=c(250,500,1000,2000,4000)


betas_clip=betas_ngd=vector("list",length=length(N))

set.seed(1020)
for(j in 1:5)
{
  betas_clip[[j]]<-matrix(nrow=ns,ncol=4)
  betas_ngd[[j]]<-matrix(nrow=ns,ncol=4)
  
  n=N[j]
  iter.j<-(50+25*(n>500)+25*(n>2000))
  for(i in 1:ns){
    x=matrix(rnorm(n*4,mean=0,sd=1),nrow=n)
    x[,1]<-rep(1,n)
    p<-invlogit(as.vector(x%*%beta_gen))
    y<-rbinom(n=n,size=1,prob=p)
    out_ngd<-NGD_logistic(x=x,y=y,private=T,maxiter=iter.j,mu=0.5*sqrt((iter.j+2)/iter.j),
                          stepsize=1,stopping=0,mnorm=sqrt(2),suppress.inference=TRUE)
    #privacy budget in the above is calibrated to make the estimation process 0.5-GDP (without allocating privacy
    #budget to inference, to be comparable to the clipping version)
    out_clip<-clipped_logistic(x=x,y=y,private=T,maxiter=iter.j,mu=0.5,c_grad=1,stepsize=1,stopping=0)

    betas_clip[[j]][i,]<-out_clip$beta
    betas_ngd[[j]][i,]<-out_ngd$beta
  }
}

#############################################
# to generate Figure [3]
#############################################

par(mfrow=c(1,2))
boxplot(betas_ngd[[1]][,4],betas_ngd[[2]][,4],betas_ngd[[3]][,4],betas_ngd[[4]][,4],betas_ngd[[5]][,4],
        names=N, xlab="sample size",ylab=expression(paste("Estimate of ", beta[4]," parameter")), ylim=c(-0.95,2.05),
        main="Logistic Regression, Mallows Weights")
abline(h=beta_gen[4],col=2)
boxplot(betas_clip[[1]][,4],betas_clip[[2]][,4],betas_clip[[3]][,4],betas_clip[[4]][,4], betas_clip[[5]][,4],
        names=N, xlab="sample size",ylab=expression(paste("Estimate of ", beta[4]," parameter")), ylim=c(-0.95,2.05),
        main="Logistic Regression, Gradient Clipping")
abline(h=beta_gen[4],col=2)



