source("PrivateOptimization/src/extra/self_concordant_optimization.R")

beta_gen=c(10,-11,9,-11)
s0=2
n<-10000

set.seed(260)
x=matrix(rnorm(n*4,mean=0,sd=2),nrow=n)
x[,1]<-rep(1,n)
y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))
out_newton<-SCNewton(x=x,y=y,private=T,mu=2*sqrt(500/25),scale=F,s0=2,stepsize=0.1,maxiter=500,stopping=0)
out_ngd<-SCNGD(x=x,y=y,private=T,mu=2*sqrt(500/375),scale=F,s0=2,maxiter=500,stopping=0)


plot(log(out_ngd$gradtraj),type="l",xlim=c(0,500),ylim=c(-8,-1),xlab="Iteration",ylab="log(norm of gradient)")
lines(log(out_newton$gradtraj),col=3)
legend("topright",legend=c("Gradient Descent","Damped Newton"),lty=1,col=c(1,3),bty='n',cex=0.8)