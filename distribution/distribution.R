# fd(x,suspects="all",graph=TRUE,ties.method="min",label=FALSE,label_series)

# Arguments: x=data
#            suspects="all" by default or =c("uniform","exponential") specified alternatively.
#            graph=TRUE by default, otherwise FALSE with no graph output
#            ties.method = c("average", "first", "random", "max", "min")
#            label=FALSE by default. If label=TRUE,you have to specify the label series.

# the argument 'ties.method' determines the result at the corresponding
# indices.  The '"first"' method results in a permutation with
# increasing values at each index set of ties.  The '"random"'
# method puts these in random order whereas 
# '"average"' replaces them by their mean, and '"max"' and '"min"'
# replaces them by their maximum and minimum respectively, the
# latter being the typical "sports" ranking and the default.


rm(list=ls())
per=0.05 # control the probability for the tail
prob=0.95 # control the quantile for labelling data. We label data when its deviation is larger than the quantile.

#################################################################################################################

library(fBasics)
library(evd)
library(stats)
library(zoo)
library(lmtest)

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

fd<-function(x,label=FALSE,label_series,suspects="all",graph=TRUE,ties.method="min")
{
 if (suspects[1]=="all")
 {suspects=c("uniform","exponential","gumbel","normal","logistic","power","pareto","frechet","weibull","log-normal","fisk")}

# removes the missing values 
  x=x[!is.na(x)]
  y=label_series[!is.na(label_series)]

 for(i in 1:length(x))
 { if (x[i]<=0)

   {
    for (i in 1:length(suspects))
    {if (suspects[i]=="power"|suspects[i]=="pareto"|suspects[i]=="frechet"|suspects[i]=="weibull"|suspects[i]=="log-normal"|suspects[i]=="fisk")
    {print(paste("Nonpositive values exist.",suspects[i],"distribution incompatible"))
     suspects[i]=NA}
    }
    
   suspects=suspects[!is.na(suspects)]
  
   break

  }
}

print("                                                                                   ")
print(paste("the skewness of data:",format(skewness(x),digits=3)))
print("                                                                                   ")
 
 combine=cbind(x,y)
 x=combine[order(combine[,1],decreasing=FALSE),1]
 x=as.numeric(x)
 y=combine[order(combine[,1],decreasing=FALSE),2]
 order<-rank(x,ties.method=ties.method)
 Fhat<-(order-0.3)/(length(order)+0.4)
 distribution=array(1,dim=c(length(suspects),11))

 for(i in 1:length(suspects))

{ if (suspects[i]=="uniform")
  {newx=x
  newFhat=Fhat
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  ltail=e[1:(per*length(x))]
  rtail=e[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(e))
  
  
  

  
  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  alpha=-fit[1]/fit[2]
  beta=(1-fit[1])/fit[2]
  ks.stat=ks.test(x,'punif',alpha,beta)$statistic
 
  
  if (graph==TRUE)
  {pdf("uniform.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}
  title(main=paste("uniform distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="uniform"
  distribution[i,2]="alpha"
  distribution[i,3]=format(alpha,digits=3)
  distribution[i,4]="beta"
  distribution[i,5]=format(beta,digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}

   
  if (suspects[i]=="exponential")
  {newx=x
  newFhat=log(1-Fhat)
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(predict)-exp(newFhat)
  ltail=newe[1:(per*length(x))]
  rtail=newe[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(newe))






  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  location=fit[1]
  scale=-fit[2]
  ks.stat=ks.test(x,'pnexp',location,scale)$statistic

  
 
  if (graph==TRUE)
  {pdf("exponential.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}


  title(main=paste("exponential distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="exponential"
  distribution[i,2]="location"
  distribution[i,3]=format(fit[1],digits=2)
  distribution[i,4]="scale"
  distribution[i,5]=format(-fit[2],digits=2)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}



  if (suspects[i]=="gumbel")
  {newx=x
  newFhat=log(-log(Fhat))
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(-exp(newFhat))-exp(-exp(predict))
  ltail=newe[1:(per*length(x))]
  rtail=newe[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(newe))






  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  location=-fit[1]/fit[2]
  scale=-fit[2]
  ks.stat=ks.test(x,'pgumbel',location,scale)$statistic
 


  if (graph==TRUE)
  {pdf("gumbel.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}


  title(main=paste("gumbel distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="gumbel"
  distribution[i,2]="location"
  distribution[i,3]=format(-fit[1]/fit[2],digits=3)
  distribution[i,4]="scale"
  distribution[i,5]=format(-fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}




  if (suspects[i]=="normal")
  {newx=x
  newFhat=qnorm(Fhat) 
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  predict=lm(newFhat~newx)$fitted.values
  newe=pnorm(newFhat)-pnorm(predict)
  ltail=newe[1:(per*length(x))]
  rtail=newe[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(newe))





  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  mu=-fit[1]/fit[2]
  sigma=1/fit[2]
  ks.stat=ks.test(x,'pnorm',mu,sigma)$statistic



  if (graph==TRUE)
  {pdf("normal.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}


  title(main=paste("normal distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="normal"
  distribution[i,2]="mu"
  distribution[i,3]=format(-fit[1]/fit[2],digits=3)
  distribution[i,4]="sigma"
  distribution[i,5]=format(1/fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}




 
  if (suspects[i]=="logistic")
  {newx=x
  newFhat=log(Fhat/(1-Fhat)) 
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(newFhat)/(1+exp(newFhat))-exp(predict)/(1+exp(predict))
  ltail=newe[1:(per*length(x))]
  rtail=newe[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(newe))





  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  location=-fit[1]/fit[2]
  scale=1/fit[2]
  ks.stat=ks.test(x,'plogis',location,scale)$statistic
 


  
  if (graph==TRUE)
  {pdf("logistic.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}


  title(main=paste("logistic distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="logistic"
  distribution[i,2]="location"
  distribution[i,3]=format(-fit[1]/fit[2],digits=3)
  distribution[i,4]="scale"
  distribution[i,5]=format(1/fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}





  if (suspects[i]=="power")
  {newx=log(x)
  newFhat=log(Fhat)
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(newFhat)-exp(predict)
  ltail=newe[1:(per*length(x))]
  rtail=newe[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(newe))





  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  scale=exp(-fit[1]/fit[2])
  shape=fit[2]
  ks.stat=ks.test(x,'ppower',scale,shape)$statistic



  
  if (graph==TRUE)
  {pdf("power.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}


  title(main=paste("power distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="power"
  distribution[i,2]="scale"
  distribution[i,3]=format(exp(-fit[1]/fit[2]),digits=3)
  distribution[i,4]="shape"
  distribution[i,5]=format(fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}





  if (suspects[i]=="pareto")
  {newx=log(x)
  newFhat=log(1-Fhat)  
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(predict)-exp(newFhat)
  ltail=newe[1:(per*length(x))]
  rtail=newe[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(newe))






  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  scale=exp(-fit[1]/fit[2])
  shape=-fit[2]
  ks.stat=ks.test(x,'ppareto',scale,shape)$statistic




  if (graph==TRUE)
  {pdf("pareto.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}


  title(main=paste("pareto distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="pareto"
  distribution[i,2]="scale"
  distribution[i,3]=format(exp(-fit[1]/fit[2]),digits=3)
  distribution[i,4]="shape"
  distribution[i,5]=format(-fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}




  if (suspects[i]=="frechet")
  {newx=log(x)
  newFhat=log(-log(Fhat)) 
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(-exp(newFhat))-exp(-exp(predict))
  ltail=newe[1:(per*length(x))]
  rtail=newe[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(newe))





  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  scale=exp(-fit[1]/fit[2])
  shape=-fit[2]
  ks.stat=ks.test(x,'pfrechet',0,scale,shape)$statistic


  
  if (graph==TRUE)
  {pdf("frechet.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}


  title(main=paste("frechet distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="frechet"
  distribution[i,2]="scale"
  distribution[i,3]=format(exp(-fit[1]/fit[2]),digits=3)
  distribution[i,4]="shape"
  distribution[i,5]=format(-fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}





  if (suspects[i]=="weibull")
  {newx=log(x)
  newFhat=log(-log(1-Fhat))  
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(-exp(predict))-exp(-exp(newFhat))
  ltail=newe[1:(per*length(x))]
  rtail=newe[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(newe))




  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  scale=exp(-fit[1]/fit[2])
  shape=fit[2]
  ks.stat=ks.test(x,'pweibull',shape,scale)$statistic




  if (graph==TRUE)
  {pdf("weibull.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}


  title(main=paste("weibull distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="weibull"
  distribution[i,2]="scale"
  distribution[i,3]=format(exp(-fit[1]/fit[2]),digits=3)
  distribution[i,4]="shape"
  distribution[i,5]=format(fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}





  if (suspects[i]=="log-normal")
  {newx=log(x)
  newFhat=qnorm(Fhat) 
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  predict=lm(newFhat~newx)$fitted.values
  newe=pnorm(newFhat)-pnorm(predict)
  ltail=newe[1:(per*length(x))]
  rtail=newe[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(newe))





  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  mu=-fit[1]/fit[2]
  sigma=1/fit[2]
  ks.stat=ks.test(x,'plnorm',mu,sigma)$statistic




  if (graph==TRUE)
  {pdf("log-normal.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}


  title(main=paste("log-normal distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="log-normal"
  distribution[i,2]="mu"
  distribution[i,3]=format(-fit[1]/fit[2],digits=3)
  distribution[i,4]="sigma"
  distribution[i,5]=format(1/fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}




  
  if (suspects[i]=="fisk")
  {newx=log(x)
  newFhat=log(Fhat/(1-Fhat)) 
  fit<-coef(lm(newFhat~newx))
  e=resid(lm(newFhat~newx))
  ht=lm(log(e^2)~log(newx))
  ht.p=summary.lm(ht)$coefficients[2,4]
  dw=dwtest(lm(newFhat~newx))$statistic
  predict=lm(newFhat~newx)$fitted.values
  newe=exp(newFhat)/(1+exp(newFhat))-exp(predict)/(1+exp(predict))
  ltail=newe[1:(per*length(x))]
  rtail=newe[((1-per)*length(x)+1):length(x)]
  left=mean(abs(ltail))
  right=mean(abs(rtail))
  deviation=mean(abs(newe))




  a=format(fit[1],digits=2)
  b=format(fit[2],digits=2)
  R=format(cor(newFhat,newx),digits=3)
  scale=exp(-fit[1]/fit[2])
  shape=fit[2]
  ks.stat=ks.test(x,'pfisk',scale,shape)$statistic


  
  if (graph==TRUE)
  {
  pdf("fisk.pdf",height=5,width=7,family="Helvetica")
  plot(newx,newFhat,type="p",lwd=2,col="red",xlab="h(x)",ylab="g(Fhat)")
  abline(fit,lwd=2,col="blue",lty="dashed")
  if(label==TRUE)
  {for(j in 1:length(x))
  {if (abs(e[j])>quantile(abs(e),prob))
    text(newx[j],newFhat[j],y[j],pos=1,offset=0.4,cex=0.5)
   }}


  title(main=paste("fisk distribution:","a=",a,"","b=",b,"","correlation coeff=",R))}

  distribution[i,1]="fisk"
  distribution[i,2]="scale"
  distribution[i,3]=format(exp(-fit[1]/fit[2]),digits=3)
  distribution[i,4]="shape"
  distribution[i,5]=format(fit[2],digits=3)
  distribution[i,6]=R
  distribution[i,7]=format(ks.stat,digits=3)
  distribution[i,8]=format(left,digits=3)
  distribution[i,9]=format(right,digits=3)
  distribution[i,10]=format(deviation,digits=3)
  distribution[i,11]=abs(cor(newFhat,newx))}




 
} # "for" loop ends


mat<-distribution[order(distribution[,11],decreasing=TRUE),1:10]

dimnames(mat)[[2]]=c("distribution","parameter1","","parameter2","","rho","ks.stat","left-tail","right-tail","mean of absolute deviation")

mat



} # function ends


#################################################################################################################

library(foreign)
data=read.dta("namedist.dta")
attach(data)

x=number

y=name

y=name[1:100]

x=rnorm(100,10,1)


#sink("result.txt")

fd(x,label=TRUE,label_series=y)

#sink()

graphics.off()



  