source("src/extra/NGD_linear_regression_extra.R")


library(MASS) #for drawing samples from a multivariate normal distribution
library(ellipse) #for plotting elliptical confidence regions

M<-c(1,0.5,0,0.5,2,0,0,0,2)
M<-matrix(M,nrow=3,byrow=T) #covariance matrix for predictors

means<-c(0,0,0)

n<-1000 #number of data points

ns<-1000 #number of repetitions

beta_gen<-c(1,1,1,1)
p<-length(beta_gen)
s0<-2


##############################
########  Outline  ###########
##############################
# 1. Simulate data (X,y) from the linear model y=XB+e
# 2. Privately estimate beta and its covariance from the simulated data, with and without the small-sample correction
#    to the covariance matrix
# 3. Repeat steps 1-2 many times
# 4. The estimated covariance matrix for (beta_2,beta_3), together with the estimated (beta_2,beta_3) values,
#    specify an elliptical confidence region for the pair (beta_2,beta_3). Calculate the "average" confidence
#    region by averaging the covariance matrices obtained across all repetitions.
# 5. Plot all pairs of estimates (beta_2,beta_3) obtained across all repetitions.
# 6. Superimpose the "average" confidence region (one using the corrected covariance, one using the non-corrected)


### Note: the plot generated below illustrates performance in the case of a non-private stopping rule
##############################

covmat23=covmat23NC=matrix(0,nrow=2,ncol=2)
center23<-rep(0,2)
betamat<-matrix(NA,nrow=ns,ncol=p)
quadform23=quadform23NC=rep(NA,ns)

set.seed(3325)

for(i in 1:ns){
  x<-mvrnorm(n=n,mu=means,Sigma=M) #generate predictor values with specified covariance matrix
  x<-cbind(rep(1,n),x)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0)) #generate response values according to linear model
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=55,mnorm=sqrt(2),stopping="non-private",suppress.inference=FALSE)
  
  betamat[i,]<-out$beta
  
  covmat23<-covmat23+(out$variances[2:3,2:3])/ns
  covmat23NC<-covmat23NC+(out$sandwich[2:3,2:3]/n)/ns
  center23<-center23+(out$beta[2:3])/ns
  
  #utility for checking whether the true (beta_2,beta_3) lies in the estimated confidence region
  #v23<-out$beta[2:3]-beta_gen[2:3]
  #quadform23[i]<-t(v23)%*%solve(out$variances[2:3,2:3])%*%v23
  #quadform23NC[i]<-t(v23)%*%solve(out$sandwich[2:3,2:3]/n)%*%v23

  
}

# utility for assessing coverage- i.e. the fraction of trials for which
# the estimated confidence region covered the true  (beta_2,beta_3)
#t<-(-2)*log(0.05)
#coverage23<-sum(quadform23<t)/ns
#coverage23NC<-sum(quadform23NC<t)/ns

ellipse23<-ellipse(covmat23,centre=center23)
ellipse23NC<-ellipse(covmat23NC,centre=center23)


##########################################
### THIS IS NOT THE PLOT IN THE PAPER#####
##########################################
par(mfrow=c(1,1))
plot(betamat[,2],betamat[,3],xlim=c(0.7,1.3),ylim=c(0.7,1.25),xlab=expression(paste(beta)[2]),ylab=expression(paste(beta)[3]), main=expression(paste("1-GDP Estimates of (",beta[2],",",beta[3],")")),pch="*")
lines(ellipse23[,1],ellipse23[,2],col=2)
lines(ellipse23NC[,1],ellipse23NC[,2],col=2,lty=3)
points(x=1,y=1,col=3,pch=16)
legend("bottomleft",legend=c("True parameter value","Avg 95% Conf Region, corrected","Avg 95% Conf Region, uncorrected"),lty=c(NULL,3,1),col=c(3,2,2),pch=c(16,NULL,NULL),bty='n',cex=0.8)






################################################################
#### this version uses a PRIVATE stopping rule, unlike the above
####    THIS ***IS*** THE PLOT IN THE PAPER      ###############
################################################################

covmat23=covmat23NC=matrix(0,nrow=2,ncol=2)
center23<-rep(0,2)
betamat<-matrix(NA,nrow=ns,ncol=p)
quadform23=quadform23NC=rep(NA,ns)

set.seed(3325)

for(i in 1:ns){
  x<-mvrnorm(n=n,mu=means,Sigma=M) #generate predictor values with specified covariance matrix
  x<-cbind(rep(1,n),x)
  y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0)) #generate response values according to linear model
  out=NGD.Huber(x=x,y=y,private=T,scale=T,mu=1,maxiter=55,mnorm=sqrt(2),stopping="private",suppress.inference=FALSE)
  
  betamat[i,]<-out$beta
  
  covmat23<-covmat23+(out$variances[2:3,2:3])/ns
  covmat23NC<-covmat23NC+(out$sandwich[2:3,2:3]/n)/ns
  center23<-center23+(out$beta[2:3])/ns
  
  #utility for checking whether the true (beta_2,beta_3) lies in the estimated confidence region
  #v23<-out$beta[2:3]-beta_gen[2:3]
  #quadform23[i]<-t(v23)%*%solve(out$variances[2:3,2:3])%*%v23
  #quadform23NC[i]<-t(v23)%*%solve(out$sandwich[2:3,2:3]/n)%*%v23
  
}

# utility for assessing coverage- i.e. the fraction of trials for which
# the estimated confidence region covered the true  (beta_2,beta_3)
#t<-(-2)*log(0.05)
#coverage23<-sum(quadform23<t)/ns
#coverage23NC<-sum(quadform23NC<t)/ns

ellipse23<-ellipse(covmat23,centre=center23)
ellipse23NC<-ellipse(covmat23NC,centre=center23)

par(mfrow=c(1,1))
plot(betamat[,2],betamat[,3],xlim=c(0.67,1.38),ylim=c(0.68,1.22),xlab=expression(paste(beta)[2]),ylab=expression(paste(beta)[3]), main=expression(paste("1-GDP Estimates of (",beta[2],",",beta[3],")")),pch="*")
lines(ellipse23[,1],ellipse23[,2],col=2)
lines(ellipse23NC[,1],ellipse23NC[,2],col=2,lty=3)
points(x=1,y=1,col=3,pch=16)
legend("bottomleft",legend=c("True parameter value","Avg 95% Conf Region, corrected","Avg 95% Conf Region, uncorrected"),lty=c(NULL,3,1),col=c(3,2,2),pch=c(16,NULL,NULL),bty='n',cex=0.8)

