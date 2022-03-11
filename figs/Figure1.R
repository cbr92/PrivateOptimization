source("NGD_linear_regression_extra.R")


beta_gen=c(1,1,1,1)
s0=2
n=1000


set.seed(9964)
x=matrix(rnorm(n*4,mean=0,sd=s0),nrow=n)
x[,1]<-rep(1,n)
y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))

outNP<-NGD.Huber(x,y,private=F,scale=T,maxiter=100,stopping=0,mnorm=sqrt(2),suppress.inference=TRUE)
NP_beta2<-outNP$beta2_traj

out_mu2<-NGD.Huber(x,y,private=T,scale=T,mu=2*sqrt(102/60),maxiter=100,stopping=0,mnorm=sqrt(2),suppress.inference=TRUE)
#mu above is calculated such that the estimation process consumes privacy budget of 2 over the first 60 iterations
private_mu2<-out_mu2$beta2_traj
out_mu05<-NGD.Huber(x,y,private=T,scale=T,mu=0.5*sqrt(102/60),maxiter=100,stopping=0,mnorm=sqrt(2),suppress.inference=TRUE)
private_mu05<-out_mu05$beta2_traj


par(mfrow=c(1,2))

plot(NP_beta2,ylab="Parameter Estimate",xlab="Iteration",ylim=c(0,1.3),type="l",main=expression(paste("Estimate of ",beta)[2]))
lines(private_mu05,col=2)
lines(private_mu2,col=4)
abline(v=60,col="gray",lty=2)
legend("bottomright",legend=c(expression(paste("private, ",mu,"=0.5")),expression(paste("private, ",mu,"=2")),"non-private"),lty=1,col=c(2,4,1),bty="n",cex=0.8)


plot(log(outNP$nonpriv_gradtraj),type="l",xlab="Iteration",ylab="log(norm of gradient)",main=expression(paste("Gradient Estimate Trajectories")))
lines(log(out_mu05$nonpriv_gradtraj),col=2)
lines(log(out_mu2$nonpriv_gradtraj),col=4)
abline(v=60,col="gray",lty=2)
legend("bottom",legend=c(expression(paste("private, ",mu,"=0.5")),expression(paste("private, ",mu,"=2")),"non-private"),lty=1,col=c(2,4,1),bty="n",cex=0.8)



