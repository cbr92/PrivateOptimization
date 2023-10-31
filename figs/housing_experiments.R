#############################################
########## CONTENTS #########################
#############################################

## 1. preprocessing of housing data set
## 2. using NGD to fit linear regression model to full data set
##    and generate Table 1
## 3. calculating average p values across repeated subsamples
##    and generage Figure [12]
## 4. calculate parameter error across repeated subsamples
##    and generate Figure [13]


#############################################
##########    1.    #########################
#############################################

#############################################
##    Title: housing.R (code), housing.Rda (data set)
##    Author: Jing Lei
##    Date: 2011
##    Code version: updated 2011, Apr 02
##    Availability: <https://www.dropbox.com/s/aakaj0m32vkiqzf/Archive.zip?dl=0>;
##                  <http://www.stat.cmu.edu/~jinglei/code.shtml>
#############################################
#load("~/Dropbox/Privacy:Casey/R code/code_Lei/housing.Rda") #data set

#housing.raw <- housing
#remove the variables we aren't going to use
#remove.ls <- c("city", "zip", "street", "long", "lat", "quality", "match")
#for (name in remove.ls) {
#  housing[, name] <- NULL
#}
#housing <- subset(housing.raw, !is.na(price)) #4 NA's removed
#housing <- subset(housing, !is.na(year)) #60879 NA's removed
#housing$date <- as.numeric(substr(housing$date, 1, 4)) #extract just the year from the date variable
#housing$date <- housing$date - 2003

#housing$county <- gsub(" County", "", as.factor(as.character(housing$county)))

#regroup the counties
#housing$county.cmb <- housing$county
#housing$county.cmb[housing$county.cmb %in% c("Marin", "San Francisco", "San Mateo")] <- "MSS"
#housing$county.cmb[housing$county.cmb %in% c("Napa", "Sonoma")] <- "NS"

#housing <- subset(housing, !is.na(bsqft)) #777 NA's removed

#housing <- housing[, c("price", "bsqft", "date", "county.cmb")] #keep only the variables we ultimately want
#housing$county.cmb<-as.factor(housing$county.cmb) #convert the county variable from character to factor
#############################################
#########  end of cited material  ###########
#############################################

#housing_mat<-model.matrix(price~.,data=housing) #convert the county variable to one-hot encoding
#colnames(housing_mat)<-c("Intercept","bsqft","date","Contra Costa", "MSS", "NS", "Santa Clara", "Solano")

#housing<-data.frame(priv.data$price,housing_mat[,2:8])
#save(housing,file="housing_data_2.Rdata")

load("~/Documents/privacy/housing_data_2.Rdata")

#### center numeric variables around their median and divide by their median absolute deviation

scaled.housingmat<-housing
for(i in 1:3){ 
  scaled.housingmat[,i]<-(housing[,i]-median(housing[,i]))/mad(housing[,i])
}

xmat<-as.matrix(scaled.housingmat[,2:8])
xmat<-cbind(rep(1,length=dim(xmat)[1]),xmat)
yvec<-scaled.housingmat[,1]

rm(housing,scaled.housingmat)



#############################################
##########    2.    #########################
#############################################

source("src/NGD_linear_regression.R")
set.seed(929)
housereg<-NGD.Huber(x=xmat,y=yvec,private=T,scale=T,mu=0.25,maxiter=100,beta0=rep(0,8),stopping=0,mnorm=sqrt(2))

coefficients<-housereg$beta
std_error<-sqrt(diag(housereg$variances))
z_value<-coefficients/std_error
pvalues<-2*pnorm(abs(z_value),lower.tail=FALSE)
summary.mat<-cbind(coefficients,std_error,z_value,pvalues)
colnames(summary.mat)<-c("Value","Std. Error", "z value", "p value")

summary.mat #table 1




#############################################
##########    3.    #########################
#############################################


N<-c(2000,4000,7500,10000,20000)
ns<-200

n<-dim(xmat)[1]
p<-dim(xmat)[2]

pvals<-vector("list",length=length(N))

set.seed(1109)

for(j in 1:length(N)){
  
  pvals[[j]]<-matrix(nrow=ns,ncol=p)
  
  for(i in 1:ns){
    idx<-sample(1:n,size=N[j],replace=TRUE)
    x<-xmat[idx,]
    y<-yvec[idx]
    
    out<-NGD.Huber(x=x,y=y,private=T,scale=T,mu=0.25,maxiter=100,beta0=rep(0,p),stopping=0,mnorm=sqrt(2))
    zscores<-out$beta/sqrt(diag(out$variances))
    pvals[[j]][i,]<-2*pnorm(abs(zscores),lower.tail=FALSE)

  }
}


avg_pvals<-matrix(nrow=length(N),ncol=p)
avg_pvals[1,]<-colMeans(pvals[[1]])
avg_pvals[2,]<-colMeans(pvals[[2]])
avg_pvals[3,]<-colMeans(pvals[[3]])
avg_pvals[4,]<-colMeans(pvals[[4]])
avg_pvals[5,]<-colMeans(pvals[[5]])

#save(N,ns,avg_pvals,pvals,file="housing_pvalues_new.Rdata")

#### one plot with all the variables

par(mfrow=c(1,2))

plot(1:5,avg_pvals[,1],col=1,xaxt='n',xlab="sample size",ylab="average p-value",main="",ylim=c(0,0.4))
axis(1,at=1:5,labels=N)
points(1:5,avg_pvals[,2],col="red1",pch=3)
points(1:5,avg_pvals[,3],col="purple1",pch=3)
points(1:5,avg_pvals[,4],col="grey50",pch=2)
points(1:5,avg_pvals[,5],col="springgreen2",pch=2)
points(1:5,avg_pvals[,6],col="steelblue2",pch=2)
points(1:5,avg_pvals[,7],col="tomato1",pch=2)
points(1:5,avg_pvals[,8],col="hotpink",pch=2)
abline(h=0.05,col="gray")
legend("topright",pch=c(1,2,3),col=c(1,1,1),legend=c("Intercept","County indicators","Numerical variables"),
       cex=0.8)

plot(1:5,log(avg_pvals[,1]),col=1,xaxt='n',xlab="sample size",ylab="log(average p-value)",main="",ylim=c(-10,0))
axis(1,at=1:5,labels=N)
points(1:5,log(avg_pvals[,2]),col="red1",pch=3)
points(1:5,log(avg_pvals[,3]),col="purple1",pch=3)
points(1:5,log(avg_pvals[,4]),col="grey50",pch=2)
points(1:5,log(avg_pvals[,5]),col="springgreen2",pch=2)
points(1:5,log(avg_pvals[,6]),col="steelblue2",pch=2)
points(1:5,log(avg_pvals[,7]),col="tomato1",pch=2)
points(1:5,log(avg_pvals[,8]),col="hotpink",pch=2)
abline(h=log(0.05),col="gray")
legend("bottomleft",pch=c(1,2,3),col=c(1,1,1),legend=c("Intercept","County indicators","Numerical variables"),
       cex=0.8)









#############################################
##########    4.    #########################
#############################################

n<-dim(xmat)[1]
p<-dim(xmat)[2]

BETA<-NGD.Huber(x=xmat,y=yvec,private=F,scale=T,maxiter=10000,beta0=rep(0,p))
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
    
    out<-NGD.Huber(x=x,y=y,private=T,scale=T,mu=0.25,maxiter=(50+25*(N[j]>500)+25*(N[j]>2000)),beta0=rep(0,p),stopping=0,suppress.inference = TRUE)
    out.np<-NGD.Huber(x=x,y=y,private=F,scale=T,maxiter=500,beta0=rep(0,p),stopping=0,suppress.inference=TRUE) #was originally used for housing plots
    
    betas_private[i,j]<-sqrt(sum((out$beta-BETA)^2))
    betas_np[i,j]<-sqrt(sum((out.np$beta-BETA)^2))
  }
}

### this is the plot for figure [13] the paper
par(mfrow=c(1,1))
plot(1:length(N),colMeans(betas_private),xlab="sample size",ylab="Error",xaxt='n',
     ylim=c(0,1),pch=16,col="cornflowerblue")
axis(1,at=1:length(N),labels=N)
lines(1:length(N),colMeans(betas_private),col="cornflowerblue")
points(1:length(N),colMeans(betas_np),col=1,pch=16)
lines(1:length(N),colMeans(betas_np),col=1)
#points(1:length(N),colMeans(betas_np_2),col=3)
legend("topright",bty='n',legend=c(expression(paste("private, ", mu,"=0.25")),"non-private"),col=c("cornflowerblue",1),lty=1,pch=16,cex=0.9)

#save(N,ns,BETA,betas_private,betas_np,file="param_error_new.Rdata")


### same plot but log scale on y axis (this is not in the paper)
plot(1:length(N),log(colMeans(betas_private)),xlab="sample size",ylab="Error",xaxt='n',
     ylim=c(-4,0),pch=16,col="cornflowerblue")
axis(1,at=1:length(N),labels=N)
lines(1:length(N),log(colMeans(betas_private)),col="cornflowerblue")
points(1:length(N),log(colMeans(betas_np)),col=1,pch=16)
lines(1:length(N),log(colMeans(betas_np)),col=1)
#points(1:length(N),colMeans(betas_np_2),col=3)
legend("topright",bty='n',legend=c("private","non-private"),col=c("cornflowerblue",1),lty=1,pch=16,cex=0.9)

