
source("src/extra/Newton_linear_regression_extra.R")
source("src/extra/NGD_linear_regression_extra.R")

beta_gen=c(1,1,1,1)
s0=2
n<-1000

set.seed(505)
x=matrix(rnorm(n*4,mean=0,sd=2),nrow=n)
x[,1]<-rep(1,n)
y=as.vector(x%*%beta_gen+rnorm(n,mean=0,sd=s0))

#non-private algorithms
outNP<-NoisyNewton(x=x,y=y,private=F,scale=F,s0=2,stopping=0,suppress.inference=TRUE)
out25NP<-NoisyNewton(x=x,y=y,private=F,scale=F,s0=2,stepsize=0.25,stopping=0,suppress.inference=TRUE)
gdNP<-NGD.Huber(x=x,y=y,private=F,scale=F,s0=2,stopping=0,suppress.inference=TRUE)

#private algorithms
outDP<-NoisyNewton(x=x,y=y,private=T,mu=2*sqrt(101/11),scale=F,s0=2,stopping=0,suppress.inference=TRUE)
out25DP<-NoisyNewton(x=x,y=y,private=T,mu=2*sqrt(101/21),scale=F,s0=2,stepsize=0.25,stopping=0,suppress.inference=TRUE)
gdDP<-NGD.Huber(x=x,y=y,private=T,mu=2*sqrt(101/71),scale=F,s0=2,stopping=0,suppress.inference=TRUE)


par(mfrow=c(1,2))

plot(log(outNP$nonpriv_gradtraj),type="l",xlab="Iteration",ylab="log(norm of gradient)",ylim=c(-10,-1),xlim=c(1,30)) #data goes to 100
lines(log(gdNP$nonpriv_gradtraj),col=2)
lines(log(out25NP$nonpriv_gradtraj),col=4)
lines(log(outDP$nonpriv_gradtraj),lty=2,col=1)
lines(log(gdDP$nonpriv_gradtraj),col=2,lty=2)
lines(log(out25DP$nonpriv_gradtraj),col=4,lty=2)
legend("bottomright",legend=c("Gradient Descent",expression(paste("Damped Newton, ",eta,"=0.25")),  
                              "Pure Newton","resp. private version"),lty=c(1,1,1,2),col=c(2,4,1,"gray"),cex=0.8)


plot(log(outNP$nonpriv_gradtraj),type="l",xlab="Iteration",ylab="log(norm of gradient)",ylim=c(-10,-1)) #data goes to 100
lines(log(gdNP$nonpriv_gradtraj),col=2)
lines(log(out25NP$nonpriv_gradtraj),col=4)
lines(log(outDP$nonpriv_gradtraj),lty=2,col=1)
lines(log(gdDP$nonpriv_gradtraj),col=2,lty=2)
lines(log(out25DP$nonpriv_gradtraj),col=4,lty=2)
legend("bottomright",legend=c("Gradient Descent",expression(paste("Damped Newton, ",eta,"=0.25")),  
                              "Pure Newton","resp. private version"),lty=c(1,1,1,2),col=c(2,4,1,"gray"),cex=0.8)

