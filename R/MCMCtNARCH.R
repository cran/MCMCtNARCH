#######################################################
# 28/10/2009 #
# Estimate the Box-Cox GARCH(1,1) model with t errors #
# using the random-walk Metropolis-Hastings algorithm #
#######################################################
rm(list=ls(all=TRUE))

##mcmctnarch

mctnarch=function(xdata, m)
{

# log posterior of parameters
log.post<-function(xp)
{
a1=exp(xp[1])/(1+exp(xp[1]))
b1=exp(xp[2])/(1+exp(xp[2]))
delta<-(-1+exp(xp[3]))/(1+exp(xp[3]))
nv=exp(xp[4])
if(a1+b1>=1) return(-exp(20))
if(nv<=3.0) return(-exp(20))

h<-alpha<-rep(0,n)
if(abs(delta)<=exp(-10))
{
alpha0<-log(h0)
alpha[1]=a1*log(y0*y0)+b1*alpha0
for(i in 2:n)
{
alpha[i]=a1*log(y[i-1]*y[i-1])+b1*alpha[i-1]
}
for(i in 1:n) h[i]=exp(alpha[i])
}
else
{
alpha0=(exp(delta*log(h0))-1)/delta
tem2=(exp(delta*log(y0*y0))-1)/delta
alpha[1]=a1*tem2+b1*alpha0
for(i in 2:n)
{
tem2=(exp(delta*log(y[i-1]*y[i-1]))-1)/delta
alpha[i]=a1*tem2+b1*alpha[i-1]
}
accept<-1
for(i in 1:n)
{
test=1+delta*alpha[i]
if(test<=0.0) accept=0
}
if(accept>0)
{
for(i in 1:n) h[i]=exp(log(1+delta*alpha[i])/delta)
}
if(accept==0) return(-exp(20))
}

#log-likelihood of tgarch
tem2<-dt(y/sqrt(h),df=nv)/sqrt(h)
#tem2<-dnorm(y/sqrt(h),mean=0,sd=1)/sqrt(h)
logpost<-sum(log(tem2))

#prior of nv
tem2<-dnorm(nv,mean=8,sd=9)
logpost=logpost+log(tem2)

#prior of a1: uniform(0,1)
#prior of b1: uniform(0,1)

#log Jacobi
logpost<-logpost+log(a1-a1*a1)+log(b1-b1*b1)
logpost<-logpost+log(delta+1)+log(1-delta)-log(2)
logpost<-logpost+xp[4]

return(logpost)
}
# Randon-walk Metropolis-Hastings Algorithm
rw.mh<-function(xp,lnpost)
{
fa=lnpost
rn=rnorm(4,mean=0,sd=1)
rn=rn/sqrt(sum(rn*rn))
tem2=xp+tune*rn
fb=log.post(tem2)
r=fb-fa
accept=0
if(r>0) accept=1
else
{ un=runif(1,min=0,max=1)
if(un<exp(r)) accept=1
}
if (accept==1)
{ xp=tem2
fa=fb
}
return(list(parameter=xp,accept=accept,logpost=fa))
}
#Main Program
set.seed(200)
nsize<-length(xdata)
xdata=xdata-mean(xdata)
n=length(xdata)-1
y=xdata[2:(n+1)]
y0=xdata[1] #the 1st observation is y0
h0=mean(xdata*xdata)
#h0 is the average of the squared returns.
a=0.05
b=0.9
delta=0.3
nv=10
para=numeric(4)
para[1]<-log(a/(1-a))
para[2]<-log(b/(1-b))
para[3]<-log((delta+1)/(1-delta))
para[4]<-log(nv)
lnpost<-log.post(para)

dm<-length(para)
warm<-m/10
tune<-0.15
para.matrix<-matrix(0,nr=m,nc=dm)
accept.rate<-0
for(ks in 1:warm)
{
tem<-rw.mh(para,lnpost)
para<-tem$parameter
accept.rate<-accept.rate+tem[[2]]
lnpost<-tem[[3]]
}

for(ks in 1:m)
{
tem<-rw.mh(para,lnpost)
para<-tem$parameter
accept.rate<-accept.rate+tem[[2]]
lnpost<-tem[[3]]
para.matrix[ks,1]<-exp(para[1])/(1+exp(para[1]))
para.matrix[ks,2]<-exp(para[2])/(1+exp(para[2]))
para.matrix[ks,3]<-(-1+exp(para[3]))/(1+exp(para[3]))
para.matrix[ks,4]<-exp(para[4])
}
accept.rate<-accept.rate/(warm+m)
accept.rate


#Compute SIF
dm<-length(para)
num.batch<-50
siz.batch<-m/num.batch
hat.para<-vector(mode="numeric",length=dm)
for(i in 1:dm)
{
hat.para[i]=mean(para.matrix[,i])
}
sif<-vector(mode="numeric",length=dm)
for(i in 1:dm)
{
tem<-matrix(para.matrix[,i],nc=num.batch,byrow=F)
batch.mn<-vector(mode="numeric",length=num.batch)
for(j in 1:num.batch)
{
batch.mn[j]<-mean(tem[,j])
}
var.batch<-siz.batch*sum((batch.mn-hat.para[i])^2)/(num.batch-1)
var.total<-sum((tem-hat.para[i])^2)/(m-1)
sif[i]<-var.batch/var.total
sif[i]<-round(sif[i],digit=2)
}
#Simulation Inefficient Factors
#The smaller its value is, the better the convergence performance is
sif
#estimate of omega, alpha, beta
round(hat.para,digit=4)
#Acceptance rate
accept.rate/m

cat("alpha:",hat.para[1],"sd of alpha:",sd(para.matrix[,1]),"sif of alpha:",sif[1],"\n")
cat("beta:",hat.para[2],"sd of beta:",sd(para.matrix[,2]),"sif of beta",sif[2],"\n")
cat("delta:",hat.para[3],"sd of delta:",sd(para.matrix[,3]),"sif of delta",sif[3],"\n")
cat("nv",hat.para[4],"sd of nv:",sd(para.matrix[,4]),"sif of nv",sif[4],"\n")

return(para.matrix)
}


###
plotpara<-function(para.matrix, m)
{

#Keep one draw for every 10 draws
ik=1:(m/10)
ik=ik*10
tem=para.matrix[ik,]
par(mfrow=c(4,1),mar=c(2,2,2,2))
plot(tem[,1],typ='l',lty=1,col=4,xlab="Iteration",ylab="para",main="a1")
plot(tem[,2],typ='l',lty=1,col=4,xlab="Iteration",ylab="para",main="b1")
plot(tem[,3],typ='l',lty=1,col=4,xlab="Iteration",ylab="para",main="delta")
plot(tem[,4],typ='l',lty=1,col=4,xlab="Iteration",ylab="para",main="nv")
}

###Example:
#xt<-matrix(scan(file="\\MCMCTNARCH\DATA\AUS2005.TXT"),ncol=1,byrow=T)
#xt<-matrix(scan(file="C:\\aus2005.txt"),ncol=1,byrow=T)
#m=10000
#mc=mctnarch(xt, m)
#plotpara(mc,m)

