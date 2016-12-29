# test(true.distribution,parameter1,parameter2,simulation,observation,seed,start,end)

# true.distribution=c("uniform","exponential","gumbel","normal","logistic","power","pareto","frechet","weibull","log-normal","fisk")
# parameter1,parameter2 specified randomly, but make sure the distribtution does not generate negative values 
# so that 11 distributions can be tested
# simulation times=100 by default
# observation size=100 by default
# By default, use set.seed(1) for random number generation , but seed can be replaced with another integer
# start,end are used to write latex file and set FLASE by default


rm(list=ls())

#################################################################################################################

library(fBasics)
library(evd)
library(stats)

# pareto

rpareto=function(n,scale=1,shape=1) scale*(1-runif(n))^(-1/shape)
ppareto=function(x, scale=1, shape=1)  1-(scale/x)^shape

# power

rpower=function(n,scale=1,shape=1) scale*(runif(n))^(1/shape)
ppower=function(x, scale=1, shape=1)  (x/scale)^shape

# fisk

rfisk=function(n,scale=1,shape=1) scale*(runif(n)/(1-runif(n)))^(1/shape)
pfisk=function(x, scale=1, shape=1)  1/(1+(scale/x)^shape)
# exponential with location parameter
rnexp=function(n,location=0,scale=1) location-scale^(-1)*log(1-runif(n))
pnexp=function(x, location=0, scale=1)  1-exp(-scale*(x-location))


#################################################################################################################

fd<-function(x,suspects="all",graph=FALSE)
{if (suspects[1]=="all")
 {suspects=c("uniform","exponential","gumbel","normal","logistic","power","pareto","frechet","weibull","log-normal","fisk")}

 x=sort(x)
 per=0.05
 order<-rank(x)
 Fhat<-(order-0.3)/(length(order)+0.4)
 distribution=array(1,dim=c(length(suspects),9))

 for(i in 1:length(suspects))

{ if (suspects[i]=="uniform")
  {newx=x
  newFhat=Fhat
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  deviation=mean(abs(e))
  

  
  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  alpha=-fit[1]/fit[2]
  beta=(1-fit[1])/fit[2]
  ks.stat=ks.test(x,'punif',alpha,beta)$statistic
  ks.p=ks.test(x,'punif',alpha,beta)$p.value
  
  if (graph==TRUE)
  {pdf("uniform.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("uniform distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="uniform"
  distribution[i,2]="alpha"
  distribution[i,3]=format(alpha,digits=3)
  distribution[i,4]="beta"
  distribution[i,5]=format(beta,digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}

   
  if (suspects[i]=="exponential")
  {newx=x
  newFhat=log(1-Fhat)
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(predict)-exp(newFhat)
  deviation=mean(abs(newe))






  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  location=fit[1]
  scale=-fit[2]
  ks.stat=ks.test(x,'pnexp',location,scale)$statistic
  ks.p=ks.test(x,'pnexp',location,scale)$p.value

  
 
  if (graph==TRUE)
  {pdf("exponential.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("exponential distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="exponential"
  distribution[i,2]="location"
  distribution[i,3]=format(fit[1],digits=2)
  distribution[i,4]="scale"
  distribution[i,5]=format(-fit[2],digits=2)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}



  if (suspects[i]=="gumbel")
  {newx=x
  newFhat=log(-log(Fhat))
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(-exp(newFhat))-exp(-exp(predict))
  deviation=mean(abs(newe))






  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  location=-fit[1]/fit[2]
  scale=-fit[2]
  ks.stat=ks.test(x,'pgumbel',location,scale)$statistic
  ks.p=ks.test(x,'pgumbel',location,scale)$p.value



  if (graph==TRUE)
  {pdf("gumbel.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("gumbel distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="gumbel"
  distribution[i,2]="location"
  distribution[i,3]=format(-fit[1]/fit[2],digits=3)
  distribution[i,4]="scale"
  distribution[i,5]=format(-fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}




  if (suspects[i]=="normal")
  {newx=x
  newFhat=qnorm(Fhat) 
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  predict=lm(newFhat~newx)$fitted.values
  newe=pnorm(newFhat)-pnorm(predict)
  deviation=mean(abs(newe))





  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  mu=-fit[1]/fit[2]
  sigma=1/fit[2]
  ks.stat=ks.test(x,'pnorm',mu,sigma)$statistic
  ks.p=ks.test(x,'pnorm',mu,sigma)$p.value


  if (graph==TRUE)
  {pdf("normal.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("normal distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="normal"
  distribution[i,2]="mu"
  distribution[i,3]=format(-fit[1]/fit[2],digits=3)
  distribution[i,4]="sigma"
  distribution[i,5]=format(1/fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}




 
  if (suspects[i]=="logistic")
  {newx=x
  newFhat=log(Fhat/(1-Fhat)) 
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(newFhat)/(1+exp(newFhat))-exp(predict)/(1+exp(predict))
  deviation=mean(abs(newe))





  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  location=-fit[1]/fit[2]
  scale=1/fit[2]
  ks.stat=ks.test(x,'plogis',location,scale)$statistic
  ks.p=ks.test(x,'plogis',location,scale)$p.value


  
  if (graph==TRUE)
  {pdf("logistic.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("logistic distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="logistic"
  distribution[i,2]="location"
  distribution[i,3]=format(-fit[1]/fit[2],digits=3)
  distribution[i,4]="scale"
  distribution[i,5]=format(1/fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}





  if (suspects[i]=="power")
  {newx=log(x)
  newFhat=log(Fhat)
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(newFhat)-exp(predict)
  deviation=mean(abs(newe))





  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  scale=exp(-fit[1]/fit[2])
  shape=fit[2]
  ks.stat=ks.test(x,'ppower',scale,shape)$statistic
  ks.p=ks.test(x,'ppower',scale,shape)$p.value


  
  if (graph==TRUE)
  {pdf("power.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("power distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="power"
  distribution[i,2]="scale"
  distribution[i,3]=format(exp(-fit[1]/fit[2]),digits=3)
  distribution[i,4]="shape"
  distribution[i,5]=format(fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}





  if (suspects[i]=="pareto")
  {newx=log(x)
  newFhat=log(1-Fhat)  
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(predict)-exp(newFhat)
  deviation=mean(abs(newe))






  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  scale=exp(-fit[1]/fit[2])
  shape=-fit[2]
  ks.stat=ks.test(x,'ppareto',scale,shape)$statistic
  ks.p=ks.test(x,'ppareto',scale,shape)$p.value



  if (graph==TRUE)
  {pdf("pareto.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("pareto distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="pareto"
  distribution[i,2]="scale"
  distribution[i,3]=format(exp(-fit[1]/fit[2]),digits=3)
  distribution[i,4]="shape"
  distribution[i,5]=format(-fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}




  if (suspects[i]=="frechet")
  {newx=log(x)
  newFhat=log(-log(Fhat)) 
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(-exp(newFhat))-exp(-exp(predict))
  deviation=mean(abs(newe))





  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  scale=exp(-fit[1]/fit[2])
  shape=-fit[2]
  ks.stat=ks.test(x,'pfrechet',0,scale,shape)$statistic
  ks.p=ks.test(x,'pfrechet',0,scale,shape)$p.value


  
  if (graph==TRUE)
  {pdf("frechet.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("frechet distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="frechet"
  distribution[i,2]="scale"
  distribution[i,3]=format(exp(-fit[1]/fit[2]),digits=3)
  distribution[i,4]="shape"
  distribution[i,5]=format(-fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}





  if (suspects[i]=="weibull")
  {newx=log(x)
  newFhat=log(-log(1-Fhat))  
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(-exp(predict))-exp(-exp(newFhat))
  deviation=mean(abs(newe))




  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  scale=exp(-fit[1]/fit[2])
  shape=fit[2]
  ks.stat=ks.test(x,'pweibull',shape,scale)$statistic
  ks.p=ks.test(x,'pweibull',shape,scale)$p.value



  if (graph==TRUE)
  {pdf("weibull.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("weibull distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="weibull"
  distribution[i,2]="scale"
  distribution[i,3]=format(exp(-fit[1]/fit[2]),digits=3)
  distribution[i,4]="shape"
  distribution[i,5]=format(fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}





  if (suspects[i]=="log-normal")
  {newx=log(x)
  newFhat=qnorm(Fhat) 
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  predict=lm(newFhat~newx)$fitted.values
  newe=pnorm(newFhat)-pnorm(predict)
  deviation=mean(abs(newe))





  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  mu=-fit[1]/fit[2]
  sigma=1/fit[2]
  ks.stat=ks.test(x,'plnorm',mu,sigma)$statistic
  ks.p=ks.test(x,'plnorm',mu,sigma)$p.value



  if (graph==TRUE)
  {pdf("log-normal.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("log-normal distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="log-normal"
  distribution[i,2]="mu"
  distribution[i,3]=format(-fit[1]/fit[2],digits=3)
  distribution[i,4]="sigma"
  distribution[i,5]=format(1/fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}




  
  if (suspects[i]=="fisk")
  {newx=log(x)
  newFhat=log(Fhat/(1-Fhat)) 
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(newFhat)/(1+exp(newFhat))-exp(predict)/(1+exp(predict))
  deviation=mean(abs(newe))




  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  scale=exp(-fit[1]/fit[2])
  shape=fit[2]
  ks.stat=ks.test(x,'pfisk',scale,shape)$statistic
  ks.p=ks.test(x,'pfisk',scale,shape)$p.value


  
  if (graph==TRUE)
  {pdf("fisk.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  title(main=paste("fisk distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="fisk"
  distribution[i,2]="scale"
  distribution[i,3]=format(exp(-fit[1]/fit[2]),digits=3)
  distribution[i,4]="shape"
  distribution[i,5]=format(fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(deviation,digits=3)
  distribution[i,9]=abs(cor(newFhat,newx))}




 
} # "for" loop ends




mat1<-distribution[order(distribution[,9],decreasing=TRUE),1]
mat2<-distribution[order(distribution[,8],decreasing=FALSE),1]


rho<-distribution[,9]

combine<-cbind(mat1,rho,mat2)

combine



} # function ends



#################################################################################################################
# test function

test=function(true.distribution,parameter1,parameter2,simulation=100,observation=100,seed=1,start=FALSE,end=FALSE)

{

set.seed(seed)

uniform=exponential=gumbel=normal=logistic=power=pareto=frechet=weibull=lognormal=fisk=0

n1=n2=n3=n4=n5=n6=n7=n8=n9=n10=n11=0

for (i in 1:simulation)

{ if (true.distribution=="pareto")

  {y=rpareto(observation,parameter1,parameter2)}

else


 { if (true.distribution=="power")

  {
   y=rpower(observation,parameter1,parameter2)}

else

{
  if (true.distribution=="gumbel")

  { 
  y=rgumbel(observation,parameter1,parameter2)}


else

{

  if (true.distribution=="frechet")

  {
  y=rfrechet(observation,0,parameter1,parameter2)}


else

{
  if (true.distribution=="fisk")

  {  y=rfisk(observation,parameter1,parameter2)}
 

else

{

  if (true.distribution=="normal")

  {y=rnorm(observation,parameter1,parameter2)} 

else

{
 
  if (true.distribution=="logistic")

  {y=rlogis(observation,parameter1,parameter2)} 

else

{
  
  if (true.distribution=="weibull")

  {y=rweibull(observation,parameter2,parameter1)} 

else

{
 
  if (true.distribution=="log-normal")

  {y=rlnorm(observation,parameter1,parameter2)} 

else

{
 
  if (true.distribution=="uniform")

  {y=runif(observation,parameter1,parameter2)} 

else

{  
  if (true.distribution=="exponential")

  {y=rnexp(observation,parameter1,parameter2)} 

}}}}}}}}}}

  
  xx=fd(y)
 
  uniform=as.numeric(xx[1,2])+uniform
  
  exponential=as.numeric(xx[2,2])+exponential

  gumbel=as.numeric(xx[3,2])+gumbel

  normal=as.numeric(xx[4,2])+normal

 logistic=as.numeric(xx[5,2])+logistic

 power=as.numeric(xx[6,2])+power

  pareto=as.numeric(xx[7,2])+pareto

 frechet=as.numeric(xx[8,2])+frechet

  weibull=as.numeric(xx[9,2])+weibull

  lognormal=as.numeric(xx[10,2])+lognormal

  fisk=as.numeric(xx[11,2])+fisk

#####################################################################
# This part controls the test direction!


smallest_deviation=xx[1,3]
highest_correlation=xx[1,1] 
#####################################################################
  if ((highest_correlation=="uniform") &(smallest_deviation=="uniform"))
  
  {n1=n1+1}

else{
  
  if ((highest_correlation=="exponential") &(smallest_deviation=="exponential"))
 
  {n2=n2+1}

else{

  if ((highest_correlation=="gumbel")  &(smallest_deviation=="gumbel"))
  
  {n3=n3+1}

else{

  if ((highest_correlation=="normal")  &(smallest_deviation=="normal"))

  
  {n4=n4+1}

else{

  if ((highest_correlation=="logistic") &(smallest_deviation=="logistic"))

  
  {n5=n5+1}
else{

  if ((highest_correlation=="power")   &(smallest_deviation=="power"))

 
  {n6=n6+1}

else{

  if ((highest_correlation=="pareto")  &(smallest_deviation=="pareto"))
 
  {n7=n7+1}

else{

  if ((highest_correlation=="frechet")  &(smallest_deviation=="frechet"))

  {n8=n8+1}

else{

  if ((highest_correlation=="weibull")  &(smallest_deviation=="weibull"))
  
  {n9=n9+1}
else{

  if ((highest_correlation=="log-normal")  &(smallest_deviation=="log-normal"))
  
  {n10=n10+1}
else{

  if ((highest_correlation=="fisk")  &(smallest_deviation=="fisk"))
  
  {n11=n11+1}

}}}}}}}}}}

}

suspects=c("","unif.","exp.","gumbel","normal","logistic","power","pareto","frechet","weibull","log-norm",paste("fisk","\\\\"))

sum=(n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11)/100

c0=true.distribution

c1=paste(format(uniform/simulation,digits=2),format(n1/sum,digits=2))

c2=paste(format(exponential/simulation,digits=2),format(n2/sum,digits=2))

c3=paste(format(gumbel/simulation,digits=2),format(n3/sum,digits=2))

c4=paste(format(normal/simulation,digits=2),format(n4/sum,digits=2))

c5=paste(format(logistic/simulation,digits=2),format(n5/sum,digits=2))

c6=paste(format(power/simulation,digits=2),format(n6/sum,digits=2))

c7=paste(format(pareto/simulation,digits=2),format(n7/sum,digits=2))

c8=paste(format(frechet/simulation,digits=2),format(n8/sum,digits=2))

c9=paste(format(weibull/simulation,digits=2),format(n9/sum,digits=2))

c10=paste(format(lognormal/simulation,digits=2),format(n10/sum,digits=2))

c11=paste(format(fisk/simulation,digits=2),format(n11/sum,digits=2),"\\\\")




mat=c(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)

ea="\\documentclass[11pt]{article}"

eb="\\begin{document}"

ec="\\begin{tabular}{|l|c|c|c|c|c|c|c|c|c|c|c|}"
ed="\\hline"
ee="\\end{tabular}"
ef="\\end{document}"


if (start==TRUE)
{write(ea,file="result.tex")
 write(eb,file="result.tex",append=TRUE)
 write("\\oddsidemargin  -1.21in",file="result.tex",append=TRUE)
 write(ec,file="result.tex",append=TRUE)
 write(ed,file="result.tex",append=TRUE)
 write(suspects,ncolumns=12,file="result.tex",append=TRUE,sep="&")}

write(ed,file="result.tex",append=TRUE)
write(mat,ncolumns=12,file="result.tex",append=TRUE,sep="&")

if (end==TRUE)
{write(ed,file="result.tex",append=TRUE)
 write(ee,file="result.tex",append=TRUE)
 write(ef,file="result.tex",append=TRUE)}


} # function ends

########################################################################################################################

test(true.distribution="uniform",parameter1=0,parameter2=1,simulation=100,observation=100,seed=1,start=TRUE)

test(true.distribution="exponential",parameter1=0,parameter2=2,simulation=100,observation=100,seed=1)

test(true.distribution="gumbel",parameter1=10,parameter2=2,simulation=100,observation=100,seed=1)

test(true.distribution="normal",parameter1=8,parameter2=2,simulation=100,observation=100,seed=1)

test(true.distribution="logistic",parameter1=10,parameter2=1,simulation=100,observation=100,seed=1)

test(true.distribution="power",parameter1=2,parameter2=2,simulation=100,observation=100,seed=1)

test(true.distribution="pareto",parameter1=2,parameter2=2,simulation=100,observation=100,seed=1)

test(true.distribution="frechet",parameter1=2,parameter2=2,simulation=100,observation=100,seed=1)

test(true.distribution="weibull",parameter1=3,parameter2=3,simulation=100,observation=100,seed=1)

test(true.distribution="log-normal",parameter1=0,parameter2=1,simulation=100,observation=100,seed=1)

test(true.distribution="fisk",parameter1=3,parameter2=2,simulation=100,observation=100,seed=1,end=TRUE)

  
  