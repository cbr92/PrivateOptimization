#############################################
########## CONTENTS #########################
#############################################

## 1. preprocessing of bank data set
## 2. calculate parameter error across repeated subsamples
##    and generate Figure [##]


#############################################
##########    1.    #########################
#############################################

#   Data set available via the UCI Machine Learning Repository
#       Dua, D. and Graff, C. (2019). UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. Irvine, CA: University of California, School of Information and Computer Science.
#   Version:  "bank-full.csv"
#   Accessible at:  https://archive.ics.uci.edu/ml/datasets/bank+marketing
#   Data set originally contributed by Moro, Cortez, and Rita
#       [Moro et al., 2014] S. Moro, P. Cortez and P. Rita. A Data-Driven Approach to Predict the Success of Bank Telemarketing. Decision Support Systems, Elsevier, 62:22-31, June 2014



bank<-read.csv("/Users/caseybradshaw/Downloads/bank/bank-full.csv")

bank<-bank[,names(bank) != "pdays"] #this  variable is technically undefined for over  80% of entries

xmat<-model.matrix(y~.,data=bank) #convert categorical predictors to one-hot encoding


# standardize numeric predictors

xmat[,"age"]<-(xmat[,"age"]-median(xmat[,"age"]))/mad(xmat[,"age"])
xmat[,"balance"]<-(xmat[,"balance"]-median(xmat[,"balance"]))/mad(xmat[,"balance"])
xmat[,"day"]<-(xmat[,"day"]-median(xmat[,"day"]))/mad(xmat[,"day"])
xmat[,"duration"]<-(xmat[,"duration"]-median(xmat[,"duration"]))/mad(xmat[,"duration"])
xmat[,"campaign"]<-(xmat[,"campaign"]-median(xmat[,"campaign"]))/mad(xmat[,"campaign"])
xmat[,"previous"]<-(xmat[,"previous"]-median(xmat[,"previous"])) #MAD of this variable is zero, so do not divide by it

yvec<-rep(0,length(bank$y))
yvec[bank$y=="yes"]<-1 #recode response as 1/0 rather than "yes"/"no"

rm(bank)

  


#############################################
##########    2.    #########################
#############################################

source("/Users/caseybradshaw/Documents/privacy/private_logistic_regression.R")

n<-dim(xmat)[1]
p<-dim(xmat)[2]

BETA<-NGD_logistic(x=xmat,y=yvec,private=F,maxiter=500,beta0=rep(0,p),mnorm=5,stepsize=1)
BETA<-BETA$beta


N<-c(500,750,1000,2000,5000,10000,15000)
ns<-400



betas_private=betas_np=matrix(nrow=ns,ncol=length(N))
set.seed(1109)
for(j in 1:length(N)){
  for(i in 1:ns){
    idx<-sample(1:n,size=N[j],replace=TRUE)
    x<-xmat[idx,]
    y<-yvec[idx]
    
    out<-NGD_logistic(x=x,y=y,private=T,mu=0.25,maxiter=(50+25*(N[j]>500)+25*(N[j]>2000)),beta0=rep(0,p),mnorm=5,stepsize=1,suppress.inference=TRUE)
    out.np<-NGD_logistic(x=x,y=y,private=F,maxiter=500,beta0=rep(0,p),mnorm=5,stepsize=1,suppress.inference=TRUE)
    
    betas_private[i,j]<-sqrt(sum((out$beta-BETA)^2))
    betas_np[i,j]<-sqrt(sum((out.np$beta-BETA)^2))
  }
}


plot(1:length(N),colMeans(betas_private),xlab="sample size",ylab="Error",xaxt='n',
     ylim=c(0,13),pch=16,col="cornflowerblue")
axis(1,at=1:length(N),labels=N)
lines(1:length(N),colMeans(betas_private),col="cornflowerblue")
points(1:length(N),colMeans(betas_np),col=1,pch=16)
lines(1:length(N),colMeans(betas_np),col=1)
legend("topright",bty='n',legend=c(expression(paste("private, ", mu,"=0.25")),"non-private"),col=c("cornflowerblue",1),lty=1,pch=16,cex=0.9)

#save(N,ns,BETA,betas_np,betas_private,file="bank_param_error_new.Rdata")


