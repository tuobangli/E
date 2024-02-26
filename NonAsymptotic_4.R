

#require foreach and doparallel for parallel processing
if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)

numCores <- detectCores()
#registering clusters, can set a smaller number using numCores-1 

registerDoParallel(numCores)

#require randtoolbox for random number generations
if (!require("randtoolbox")) install.packages("randtoolbox")
library(randtoolbox)
#require Rfast for faster computation
if (!require("Rfast")) install.packages("Rfast")
library(Rfast)
if (!require("gnorm")) install.packages("gnorm")
library(gnorm)
if (!require("stats")) install.packages("stats")
library(stats)
if (!require("np")) install.packages("np")
library(np)
if (!require("nloptr")) install.packages("nloptr")
library(nloptr)

factdivide<-function(n1,n2){
  decin1<-n1-floor(n1)
  if(decin1==0){decin1=1}
  decin2<-n2-floor(n2)
  if(decin2==0){decin2=1}
  n1seq<-seq(decin1,n1, by=1)
  n2seq<-seq(decin2,n2,by=1)
  all<-list(n1seq,n2seq)
  maxlen <- max(lengths(all))
  all2 <- as.data.frame(lapply(all, function(lst) c(lst, rep(1, maxlen - length(lst)))))
  division<-all2[,1]/all2[,2]
  answer<-exp(sum(log(division)))*(gamma(decin1)/gamma(decin2))
  return(answer)
}
unbiasedsd<-function (x){
  n<-length(x)
  if(n==1){
    return(1000000)
  }
  sd1<-sd(x)
  if(n==2){
    c1<-0.7978845608
  }else if (n==3){
    c1<-0.8862269255 
  }else if (n==4){
    c1<-0.9213177319 
  }else{
    c1<-sqrt(2/(n-1))*factdivide(n1=((n/2)-1),n2=(((n-1)/2)-1))
  }
  listall<-sd1*c1
  (listall)
}
correctfactor<-function (n){
  if(n==2){
    c1<-0.7978845608
  }else if (n==3){
    c1<-0.8862269255 
  }else if (n==4){
    c1<-0.9213177319 
  }else{
    c1<-sqrt(2/(n-1))*factdivide(n1=((n/2)-1),n2=(((n-1)/2)-1))
  }
  listall<-c1
  (listall)
}
#load the deterministic simulation functions of 9 common unimodal distributions
dsexp<-function (uni,scale=1) {
  sample1<-qexp(uni,rate=scale)
  sample1
}
dsRayleigh<-function (uni,scale=1) {
  sample1 <- scale * sqrt(-2 * log((uni)))
  sample1[scale <= 0] <- NaN
  rev(sample1)
}
dsnorm<-function (uni,location=0,scale=1) {
  sample1<-qnorm(uni,mean =location,sd=scale)
  sample1
}
dsLaplace<-function (uni,location=0,scale=1) {
  sample1<-location - sign(uni - 0.5) * scale * (log(2) + ifelse(uni < 0.5, log(uni), log1p(-uni)))
  sample1
}
dslogis<-function (uni,location=0,scale=1) {
  sample1<-qlogis(uni,location=location,scale=scale)
  sample1
}
dsPareto<-function (uni,shape,scale=1) {
  sample1 <- scale*((uni))^(-1/shape)
  sample1[scale <= 0] <- NaN
  sample1[shape <= 0] <- NaN
  rev(sample1)
}
dslnorm<-function (uni,location=0,scale) {
  sample1 <- qlnorm(uni,meanlog=location,sdlog = scale)
  sample1
}
dsgamma<-function (uni,shape,scale = 1) {
  sample1<-qgamma(uni,shape=shape,scale=scale)
  sample1
}
dsWeibull<-function (uni,shape, scale = 1){
  sample1<-qweibull(uni,shape=shape, scale = scale)
  sample1
}
dsgnorm<-function (uni,location,shape, scale = 1){
  sample1<-qgnorm(p=uni, mu = location, alpha = scale, beta = shape)
  sample1
}

dsbeta<-function (uni,shape1,shape2) {
  sample1<-qbeta(uni,shape1=shape1,shape2=shape2)
  sample1
}
      
dsBSSN<-function (uni,mu, sigma, nu, tau) {
  library(gamlssbssn)
  sample1<-qBSSN(uni,mu=mu,sigma=sigma,nu=nu,tau=tau,lower.tail = TRUE,log.p=FALSE)
  sample1
}

#moments for checking the accuracy of bootstrap 
moments<-function (x){
  n<-length(x)
  m1<-mean(x)
  var1<-(sum((x - m1)^2)/n)
  tm1<-(sum((x - m1)^3)/n)
  fm1<-(sum((x - m1)^4)/(n))
  listall<-c(mean=m1,variance=var1,tm=tm1,fm=fm1)
  (listall)
}
unbiasedmoments<-function (x){
  n<-length(x)
  m1<-mean(x)
  var1<-sd(x)^2
  var2<-(sum((x - m1)^2)/n)
  tm1<-(sum((x - m1)^3)/n)*(n^2/((n-1)*(n-2)))
  fm1<-(sum((x - m1)^4)/n)
  ufm1<--3*var2^2*(2*n-3)*n/((n-1)*(n-2)*(n-3))+(n^2-2*n+3)*fm1*n/((n-1)*(n-2)*(n-3))
  listall<-c(mean=m1,variance=var1,tm=tm1,fm=ufm1)
  (listall)
}

standardizedmoments<-function (x){
  n<-length(x)
  m1<-mean(x)
  var1<-sd(x)^2
  var2<-(sum((x - m1)^2)/n)
  tm1<-(sum((x - m1)^3)/n)*(n^2/((n-1)*(n-2)))
  fm1<-(sum((x - m1)^4)/n)
  ufm1<--3*var2^2*(2*n-3)*n/((n-1)*(n-2)*(n-3))+(n^2-2*n+3)*fm1*n/((n-1)*(n-2)*(n-3))
  listall<-c(mean=m1,variance=var1,skewness=tm1/((var1)^(3/2)),kurtosis=ufm1/((var1)^(2)))
  (listall)
}

biasedmoments_expected<-function (n,targetm,targetvar,targettm,targetfm){
  m1<-targetm
  var1=targetvar/(n/(n-1))
  
  tm1<-targettm/(n^2/((n-1)*(n-2)))
  
  fm1<-(targetfm+3*(targetvar/(n/(n-1)))^2*(2*n-3)*n/((n-1)*(n-2)*(n-3)))/((n^2-2*n+3)*n/((n-1)*(n-2)*(n-3)))
  
  listall<-c(mean=m1,variance=var1,tm=tm1,fm=fm1)
  (listall)
}
round_sum_preserved<-function(x,digits=0){
  x<-x*(10^digits)
  floorx<-floor(x)
  order1<-tail(order(x-floorx),round(sum(x))-sum(floorx))
  floorx[order1]<-floorx[order1]+1
  results1<-floorx/(10^digits)
  return(results1)
}
weighted_SE<-function(x,weights){
  i <- !is.na(x)
  weights <- weights[i]
  x <- x[i]
  n_eff1 <- (sum(weights))^2/(sum(weights^2))
  wSE1<-sqrt((n_eff1/(n_eff1-1) * (sum(weights*(x-weighted.mean(x,weights))^2)/sum(weights)))/n_eff1)
  return(wSE1)
}
weighted_SD<-function(x,weights){
  i <- !is.na(x)
  weights <- weights[i]
  x <- x[i]
  SD1 =sqrt(sum(weights*(x-weighted.mean(x,weights))^2)/sum(weights))
  return(SD1)
}


beta_arithematic_sequences_process<-function(function1=NULL,distype=NULL,parameters=NULL,expect1=NULL,samplesize=NULL,seed1=NULL,weight1=NULL){
  beta_U1<-0.547

  set.seed(seed1)
  dataframe_Gaussian<-c()
  
  #first calculate the expected bias based on unbiased moments
  Results1<-c(0,expect1)
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #the first distribution used to approximate the finite sample distribution is based on the arithmetic sequence
  
  length2<-samplesize
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(1,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #then, beta distribution with U-shape
  
  length2<-samplesize
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  uniran_beta_U1<-dsbeta(uni=uniran1,shape1 =beta_U1,shape2=beta_U1)
  
  x<-do.call(distype, c(list(uniran_beta_U1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(2,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  length2<-samplesize
  uniran1<-runif(length2)
  
  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(3,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  length2<-samplesize
  uniran1<-runif(length2)
  
  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(4,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)

  length2<-samplesize
  
  uniran1<-runif(length2)
  
  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(5,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)

  length2<-samplesize
  
  uniran1<-runif(length2)
  
  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(6,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)

  length2<-samplesize
  
  uniran1<-runif(length2)
  
  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(7,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)

  length2<-samplesize
  
  uniran1<-runif(length2)
  
  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(8,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)

  length2<-samplesize
  
  uniran1<-runif(length2)
  
  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(9,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)

  length2<-samplesize
  
  uniran1<-runif(length2)

  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(10,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  length2<-samplesize
  
  uniran1<-runif(length2)
  
  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(11,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  length2<-samplesize
  
  uniran1<-runif(length2)

  x<-do.call(distype, c(list(uniran1), parameters)) 
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(12,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  
  dataframe1<-rbind(dataframe_Gaussian)
  
  return(dataframe1)
}

eval_f <- function(x){
  return ((x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
             x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
           +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
           -dataframe1[1,2])^2+(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                  x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                -dataframe1[1,3])^2+(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                       x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                     +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                     -dataframe1[1,4])^2+(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                            x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                          +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                          -dataframe1[1,5])^2)
}
eval_grad_f <-function(x){
  return(c(2*dataframe_Process[1,1]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,1]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,1]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,1]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5]),
           2*dataframe_Process[1,2]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,2]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,2]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,2]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5]),
           2*dataframe_Process[1,3]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,3]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,3]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,3]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5]),
           2*dataframe_Process[1,4]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,4]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,4]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,4]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5]),
           2*dataframe_Process[1,5]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,5]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,5]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4]),
           2*dataframe_Process[1,6]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,6]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,6]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,6]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5]),
           2*dataframe_Process[1,7]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,7]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,7]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,7]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5]),
           2*dataframe_Process[1,8]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,8]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,8]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,8]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5]),
           2*dataframe_Process[1,9]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,9]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,9]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,9]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5]),
           2*dataframe_Process[1,10]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                        x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                      +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                      -dataframe1[1,2])+2*dataframe_Process[2,10]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                     x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                   +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                   -dataframe1[1,3])+2*dataframe_Process[3,10]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                                  x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                                +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                                -dataframe1[1,4])+2*dataframe_Process[4,10]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                               x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                             +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                             -dataframe1[1,5]),
           2*dataframe_Process[1,11]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                        x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                      +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                      -dataframe1[1,2])+2*dataframe_Process[2,11]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                     x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                   +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                   -dataframe1[1,3])+2*dataframe_Process[3,11]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                                  x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                                +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                                -dataframe1[1,4])+2*dataframe_Process[4,11]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                               x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                             +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                             -dataframe1[1,5]),
           2*dataframe_Process[1,12]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                        x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                      +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                      -dataframe1[1,2])+2*dataframe_Process[2,12]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                     x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                   +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                   -dataframe1[1,3])+2*dataframe_Process[3,12]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                                  x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                                +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                                -dataframe1[1,4])+2*dataframe_Process[4,12]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                               x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                             +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                             -dataframe1[1,5]))) 
}
eval_g_eq<-function(x){
  constr<-c((1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10]))^2+(x[3]-x[4])^2+(x[5]-x[6])^2+(x[7]-x[8])^2+(x[9]-x[10])^2)
  grad<-c(-2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          2*(x[3]-x[4])-2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(x[3]-x[4])-2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          2*(x[5]-x[6])-2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(x[5]-x[6])-2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          2*(x[7]-x[8])-2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(x[7]-x[8])-2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          2*(x[9]-x[10])-2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(x[9]-x[10])-2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10]))
  )
  return(list("constraints"=constr,"jacobian"=grad))
}
x0<-rep(1,12)
lb<-rep(0,12)
ub<-rep(1,12)
opts <- list("algorithm"="NLOPT_LD_SLSQP",
             "xtol_rel"=1.0e-10,"maxeval" = 1000)



simulatedbatch_Finitesample<-foreach(batchnumber =rep(5:100,20), .combine = 'rbind') %dopar% {
  library(Rfast)
  if (!require("foreach")) install.packages("foreach")
  library(foreach)
  if (!require("doParallel")) install.packages("doParallel")
  library(doParallel)
  #registering clusters, can set a smaller number using numCores-1 
  
  #require randtoolbox for random number generations
  if (!require("randtoolbox")) install.packages("randtoolbox")
  library(randtoolbox)
  #require Rfast for faster computation
  if (!require("Rfast")) install.packages("Rfast")
  library(Rfast)
  if (!require("gnorm")) install.packages("gnorm")
  library(gnorm)
  library(nloptr)
  samplesize=batchnumber
  #seed1=1
  #set.seed(batchnumber)
  distype =dsWeibull
  distypename="Weibull1"
  a=1
  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  
  step2<-1
  repeat{
    set.seed(round(runif(1,min=1,max=2)*batchnumber*step2*100))
    random_numbers <- runif(12)
    sum_random_numbers <- sum(random_numbers)
    weight1 <- random_numbers / sum_random_numbers
    random_numbers <- runif(12)
    sum_random_numbers <- sum(random_numbers)
    weight2 <- random_numbers / sum_random_numbers
    
    step1<-1
    seed11<-round(runif(1)*10000000)
    repeat{
      dataframe1<-beta_arithematic_sequences_process(function1 = moments,distype =distype,parameters=c(shape=a, scale = 1),expect1=biasedmoments_expected(n=samplesize,targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm),
                                                     samplesize=samplesize,seed1=seed11,weight1=weight1)

      dataframe_Process<-t(cbind(dataframe1[2:13,2:5]))
      
      x0<-weight2
      lb<-rep(0,12)
      ub<-rep(1,12)
      opts <- list("algorithm"="NLOPT_LD_SLSQP",
                   "xtol_rel"=1.0e-10,"maxeval" = 1000)
      
      library(nloptr)
      results1 <- nloptr( x0=x0,
                          eval_f=eval_f,
                          lb=lb,
                          ub=ub,
                          eval_grad_f=eval_grad_f,
                          eval_g_eq=eval_g_eq,
                          opts=opts)
      
      step1<-step1+1
      weight1<-results1$solution
      
      diffweights<-sum(abs(weight1-weight2))
      if (step1==100||diffweights<1e-5){
        break
      }
      weight2<-weight1
    }
    step2<-step2+1
    if ((step2>1000&&diffweights<1.0e-5&&results1$objective<1e-2)||results1$objective<1e-10){
      print(samplesize)
      break
    }
    
  }
  
  write.csv(c(samplesize=samplesize,seed1=seed11,objective=results1$objective,weights=results1$solution,distypename=distypename),file=paste("weights",distypename,batchnumber,seed11,".csv", sep = ","), row.names=FALSE)
  
  c(samplesize=samplesize,seed1=seed11,objective=results1$objective,weights=results1$solution,distypename)
}
write.csv(simulatedbatch_Finitesample,file=paste("weights_Weibull_1.csv", sep = ","), row.names=FALSE)

library(dplyr)
simulatedbatch_Finitesample<-read.csv(file=paste("weights_Weibull_1.csv", sep = ","))

simulatedbatch_Finitesample<-simulatedbatch_Finitesample[,1:15]
# Group the data by the first column and compute the mean of other columns for each group
grouped_df <- simulatedbatch_Finitesample %>%
  group_by(samplesize) %>%
  summarize_all(mean)

# View the grouped dataframe
print(grouped_df)
write.csv(grouped_df,file=paste("weights_Weibull_1_group.csv", sep = ","), row.names=FALSE)





stopCluster(numCores)
registerDoSEQ()

