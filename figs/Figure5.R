
source("PrivateOptimization/src/extra/Newton_linear_regression_extra.R")
source("PrivateOptimization/src/extra/NGD_linear_regression_extra.R")

beta_gen=c(1,1,1,1)
s0=2
n<-1000


set.seed(222)

x=matrix(rnorm(n*4,mean=0,sd=2),nrow=n)
x[,1]<-rep(1,n)
y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))

out<-NoisyNewton(x=x,y=y,private=T,mu=2*sqrt(121/9),scale=F,s0=2,maxiter=120,stopping=0,suppress.inference=TRUE)
outNP<-NoisyNewton(x=x,y=y,private=F,scale=F,s0=2,maxiter=120,stopping=0,suppress.inference=TRUE)

gd<-NGD.Huber(x=x,y=y,private=T,mu=2*sqrt(121/81),scale=F,s0=2,maxiter=120,stopping=0,suppress.inference=TRUE)
gdNP<-NGD.Huber(x=x,y=y,private=F,scale=F,s0=2,maxiter=120,stopping=0,suppress.inference=TRUE)




plot(log(out$nonpriv_gradtraj),type="l",xlim=c(1,120),ylim=c(-7,-1),ylab="log(norm of gradient)",xlab="Iteration",col=2)
points(log(gd$nonpriv_gradtraj),col=3,type="l")
lines(log(gdNP$nonpriv_gradtraj),col=4)
lines(log(outNP$nonpriv_gradtraj),col=1)
legend("topright", legend=c("Non-private Newton","Private Newton", "Non-private GD", "Private GD"),lty=1,col=c(1,2,4,3),bty='n',cex=0.8)
abline(v=8,col="gray",lty=2)
abline(v=80,col="gray",lty=2)
