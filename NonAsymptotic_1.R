

#require foreach and doparallel for parallel processing
if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)

numCores <- detectCores()
#registering clusters, can set a smaller number using numCores-1 

registerDoParallel(numCores-1)

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

beta_arithematic_sequences_process<-function(function1=NULL,expect1=NULL,samplesize=NULL,seed1=NULL,weight1=NULL){
  beta_U1<-0.547
  
  beta_n1=20.108
  beta_n2=46.761
  
  beta_ari1<-0.328
  
  beta_L1<-0.478
  beta_L2<-38.53
  
  beta_beta1<-0.369
  beta_beta2<-18.933
  
  dataframe_Gaussian<-c()
  
  #first calculate the expected bias based on unbiased moments
  Results1<-c(0,expect1)
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #the first distribution used to approximate the finite sample distribution is based on the arithmetic sequence
  
  length2<-samplesize
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  
  x<-dsnorm(uni=uniran1,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(1,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #then, beta distribution with U-shape
  
  length2<-samplesize
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  uniran_beta_U1<-dsbeta(uni=uniran1,shape1 =beta_U1,shape2=beta_U1)
  
  x<-dsnorm(uni=uniran_beta_U1,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(2,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #then, beta distribution with n-shape, left skewed
  
  length2<-samplesize
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  uniran_n1<-dsbeta(uni=uniran1,shape1 =beta_n1,shape2=beta_n2)
  
  x<-dsnorm(uni=uniran_n1,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(3,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #beta distribution with n-shape, right skewed
  
  length2<-samplesize
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  uniran_n2<-dsbeta(uni=uniran1,shape1 =beta_n2,shape2=beta_n1)
  
  x<-dsnorm(uni=uniran_n2,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(4,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #beta distribution mixed with arithematic sequence (right)
  
  length2<-samplesize
  
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  
  uniran_ari1<-dsbeta(uni=uniran1,shape1 =beta_ari1,shape2=beta_ari1)
  
  uniran_beta_ari1<-c(uniran1[1:(round(samplesize/2))],uniran_ari1[(round(samplesize/2)+1):samplesize])
  
  x<-dsnorm(uni=uniran_beta_ari1,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(5,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #beta distribution mixed with arithematic sequence (left)
  
  length2<-samplesize
  
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  
  uniran_ari1<-dsbeta(uni=uniran1,shape1 =beta_ari1,shape2=beta_ari1)
  
  uniran_beta_ari2<-c(uniran_ari1[1:(round(samplesize/2))],uniran1[(round(samplesize/2)+1):samplesize])
  
  x<-dsnorm(uni=uniran_beta_ari2,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(6,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #beta distribution right skewed
  
  length2<-samplesize
  
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  
  uniran_beta_L1<-dsbeta(uni=uniran1,shape1 =beta_L1,shape2=beta_L2)
  
  x<-dsnorm(uni=uniran_beta_L1,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(7,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #beta distribution left skewed
  
  length2<-samplesize
  
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  
  uniran_beta_L2<-dsbeta(uni=uniran1,shape1 =beta_L2,shape2=beta_L1)
  
  x<-dsnorm(uni=uniran_beta_L2,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(8,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #beta distribution left skewed and right skewed mixed
  
  length2<-samplesize
  
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  
  uniran_beta_beta1<-dsbeta(uni=uniran1,shape1 =beta_beta1,shape2=beta_beta1)
  
  uniran_beta_beta2<-dsbeta(uni=uniran1,shape1 =beta_beta2,shape2=beta_beta2)
  
  uniran_beta_betaa<-c(uniran_beta_beta1[1:(round(samplesize/2))],uniran_beta_beta2[(round(samplesize/2)+1):samplesize])
  
  x<-dsnorm(uni=uniran_beta_betaa,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(9,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #beta distribution left skewed and right skewed mixed
  
  length2<-samplesize
  
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  
  uniran_beta_beta1<-dsbeta(uni=uniran1,shape1 =beta_beta1,shape2=beta_beta1)
  
  uniran_beta_beta2<-dsbeta(uni=uniran1,shape1 =beta_beta2,shape2=beta_beta2)
  
  uniran_beta_betab<-c(uniran_beta_beta2[1:(round(samplesize/2))],uniran_beta_beta1[(round(samplesize/2)+1):samplesize])
  
  x<-dsnorm(uni=uniran_beta_betab,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(10,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  #random sequence with antivariate
  
  length2<-samplesize
  set.seed(seed1)
  ran1<-runif(length2)
  uniran<-ran1
  
  x<-dsnorm(uni=uniran,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(11,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  
  #a final complement sequence that making the overall sequences a uniform shape
  
  length2<-samplesize
  
  factor1<-(1/(weight1[1:11]))
  
  if (sum(factor1[c(factor1>1 & factor1<10000)])>0){
    weight2<-round(weight1*(max(factor1[c(factor1>1 & factor1<10000)])))
  }else{
    weight2<-rep(1,12)
  }
  
  uniranall_weighted<-c(rep(uniran1,weight2[1]),rep(uniran_beta_U1,weight2[2]),rep(uniran_n1,weight2[3]),rep(uniran_n2,weight2[4]),
                        rep(uniran_beta_ari1,weight2[5]),rep(uniran_beta_ari2,weight2[6]),rep(uniran_beta_L1,weight2[7]),rep(uniran_beta_L2,weight2[8]),rep(uniran_beta_betaa,weight2[9]),rep(uniran_beta_betab,weight2[10]),rep(uniran,weight2[11]))
  
  uniranlong<-seq(from=1/(2*length2+1), to=1-1/(2*length2+1), by=1/(2*length2+1))
  
  uniranall_Table <- as.data.frame(table(cut(uniranall_weighted,seq(0,1,1/(2*length2)))))
  
  if (weight2[12]>0){
    length3<-(length(uniranall_weighted)+length2*weight2[12])
  }else{
    length3<-length(uniranall_weighted)*((1.1))
  }
  uniranall_Table$Freq2 <- length3/((2*length2))-uniranall_Table$Freq
  
  uniranall_Table$Freq3<-uniranall_Table$Freq2/sum(uniranall_Table$Freq2)
  
  uniranall_Table$Levels<-uniranlong
  uniranall_Table$Freq4<-uniranall_Table$Freq3*length2
  
  if (min(uniranall_Table$Freq4)<0){
    uniranall_Table$Freq5 <- (uniranall_Table$Freq4)-min(uniranall_Table$Freq4)
  }else{
    uniranall_Table$Freq5<-uniranall_Table$Freq4
  }
  
  uniranall_Table$Freq6<-uniranall_Table$Freq5/(sum(uniranall_Table$Freq5)/sum(uniranall_Table$Freq4))
  uniranall_Table$Freq7<-round_sum_preserved(uniranall_Table$Freq6)
  
  uniran_complement <- rep(uniranall_Table$Levels, uniranall_Table$Freq7)
  
  x<-dsnorm(uni=uniran_complement,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(12,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  
  dataframe1<-rbind(dataframe_Gaussian)
  
  return(dataframe1)
}


beta_arithematic_random_sequences_process<-function(function1=NULL,distype=NULL,parameters=NULL,expect1=NULL,samplesize=NULL,seed1=NULL,weight1=NULL){
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




simulatedbatch_Finitesample<-foreach(batchnumber =rep(5:100,10), .combine = 'rbind') %dopar% {
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
  step2<-1
  repeat{
    weight1<-rep(runif(1),12)
    weight2<-rep(runif(1),12)
    
    step1<-1
    seed11<-round(runif(1)*10000000)
  repeat{
    
    dataframe1<-beta_arithematic_sequences_process(function1 = moments,expect1=biasedmoments_expected(n=samplesize,targetm=0,targetvar=1,targettm=0,targetfm=3),
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
    if ((step2>10000&&diffweights<1.0e-5)||results1$objective<1e-10){
      print(samplesize)
      break
    }
    
  }
  write.csv(c(samplesize=samplesize,seed1=seed11,objective=results1$objective,weights=results1$solution),file=paste("weights",batchnumber,seed11,".csv", sep = ","), row.names=FALSE)
  
  c(samplesize=samplesize,seed1=seed11,objective=results1$objective,weights=results1$solution)
}
write.csv(simulatedbatch_Finitesample,file=paste("weights_designed_Gaussian.csv", sep = ","), row.names=FALSE)
library(dplyr)
library(readr)

library(dplyr)
simulatedbatch_Finitesample<-read.csv(file=paste("weights_designed_Gaussian.csv", sep = ","))

simulatedbatch_Finitesample<-simulatedbatch_Finitesample[,1:15]
# Group the data by the first column and compute the mean of other columns for each group
grouped_df <- simulatedbatch_Finitesample %>%
  group_by(samplesize) %>%
  summarize_all(mean)

# View the grouped dataframe
print(grouped_df)
write.csv(grouped_df,file=paste("weights_designed_Gaussian_group.csv", sep = ","), row.names=FALSE)


#then use the calculated weight and seed number(corresponds to the randome sequences selected) to estimate finite sample bias
simulatedbatch_Finitesample<-as.data.frame(simulatedbatch_Finitesample)


Wmean1<-c()
Wsd1<-c()
for (samplesize2 in (5:100)){
subset_df <- simulatedbatch_Finitesample[simulatedbatch_Finitesample$samplesize %in% samplesize2, ]
all_Gaussian_sd<-c()
dataframe2<-c()
for(batch1 in (1:nrow(subset_df))){
  seed2<-as.numeric(subset_df[batch1,2])
  weight2<-as.numeric(subset_df[batch1,4:15])
  #samplesize2<-as.numeric(subset_df[batch1,1])
  dataframe2<-beta_arithematic_sequences_process(function1=sd,expect1=0,samplesize=samplesize2,seed1=seed2,weight1=weight2)
  dataframe_Process_Gaussian<-(cbind(dataframe2[2:13,2],weight2))
  all_Gaussian_sd<-rbind(all_Gaussian_sd,dataframe_Process_Gaussian)
}

Wmean0<-weighted.mean(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
Wsd0<-weighted_SD(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
Wmean1<-rbind(Wmean1,c(samplesize2,Wmean0))

Wsd1<-rbind(Wsd1,c(samplesize2,Wsd0))
}
write.csv(Wmean1,file=paste("SD_designed_Gaussian.csv", sep = ","), row.names=FALSE)
write.csv(Wsd1,file=paste("SD_sd_designed_Gaussian.csv", sep = ","), row.names=FALSE)


Wsdtrue<-c()
for (samplesize2 in (5:100)){
  factorsd<-sqrt((1-(2/(samplesize2-1))*((gamma(samplesize2/2)/gamma((samplesize2-1)/2))^2)))
  Wsdtrue<-rbind(Wsdtrue,c(samplesize2,factorsd))
  
}
write.csv(Wsdtrue,file=paste("SD_sd_true_Gaussian.csv", sep = ","), row.names=FALSE)

Wmeantrue<-c()
for (samplesize2 in (5:100)){
  
  Wmeantrue<-rbind(Wmeantrue,c(samplesize2,correctfactor(samplesize2)))
  
}
write.csv(Wmeantrue,file=paste("SD_true_Gaussian.csv", sep = ","), row.names=FALSE)

Wmeanarith<-c()
for (samplesize2 in (5:100)){
  length2<-samplesize2
  uniran1<-seq(from=1/(length2+1), to=1-1/(length2+1), by=1/(length2+1))
  
  x<-dsnorm(uni=uniran1,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)

  Wmeanarith<-rbind(Wmeanarith,c(samplesize2,sd(sortedx)))
  
}
write.csv(Wmeanarith,file=paste("SD_arith_Gaussian.csv", sep = ","), row.names=FALSE)


Wmeanmonte<-c()
Wsdmonte1<-c()
for (samplesize2 in (5:100)){
  length2<-samplesize2
  sdnorm<-c()
  for (i in (1:120)){
    x<-rnorm(length2)
    sdnorm<-c(sdnorm,sd(x))
  }
  
  Wmeanmonte<-rbind(Wmeanmonte,c(samplesize2,mean(sdnorm)))
  Wsdmonte1<-rbind(Wsdmonte1,c(samplesize2,sd(sdnorm)))
}

write.csv(Wmeanmonte,file=paste("SD_Monte_Gaussian_120.csv", sep = ","), row.names=FALSE)
write.csv(Wsdmonte1,file=paste("SD_Monte_Gaussian_120_sd.csv", sep = ","), row.names=FALSE)

simulatedbatch_Finitesample<-read.csv(file=paste("weights_beta-arith_random.csv", sep = ","))

simulatedbatch_Finitesample<-as.data.frame(simulatedbatch_Finitesample)
  
Wmean01<-c()
Wsd01<-c()
for (samplesize2 in (5:100)){
  subset_df <- simulatedbatch_Finitesample[simulatedbatch_Finitesample$samplesize %in% samplesize2, ]
  all_Gaussian_sd<-c()
  dataframe2<-c()
  for(batch1 in (1:nrow(subset_df))){
    seed2<-as.numeric(subset_df[batch1,2])
    weight2<-as.numeric(subset_df[batch1,4:15])
    #samplesize2<-as.numeric(subset_df[batch1,1])
    dataframe2<-beta_arithematic_random_sequences_process(function1=sd,distype=dsnorm,parameters=c(location=0,scale=1),expect1=0,samplesize=samplesize2,seed1=seed2,weight1=weight2)
    dataframe_Process_Gaussian<-(cbind(dataframe2[2:13,2],weight2))
    all_Gaussian_sd<-rbind(all_Gaussian_sd,dataframe_Process_Gaussian)
  }
  
  Wmean00<-weighted.mean(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
  Wsd00<-weighted_SD(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
  Wmean01<-rbind(Wmean01,c(samplesize2,Wmean00))
  
  Wsd01<-rbind(Wsd01,c(samplesize2,Wsd00))
}
write.csv(Wmean01,file=paste("SD_beta-arith_random_Gaussian_120.csv", sep = ","), row.names=FALSE)
write.csv(Wsd01,file=paste("SD_beta-arith_random_Gaussian_120_sd.csv", sep = ","), row.names=FALSE)

simulatedbatch_Finitesample1<-read.csv(file=paste("weights_beta-arith_random.csv", sep = ","))

simulatedbatch_Finitesample2<-read.csv(file=paste("weights_norm1.csv", sep = ","))

simulatedbatch_Finitesample<-as.data.frame(rbind(simulatedbatch_Finitesample1,simulatedbatch_Finitesample2))

Wmean02<-c()
Wsd02<-c()
for (samplesize2 in (5:100)){
  subset_df <- simulatedbatch_Finitesample[simulatedbatch_Finitesample$samplesize %in% samplesize2, ]
  all_Gaussian_sd<-c()
  dataframe2<-c()
  for(batch1 in (1:nrow(subset_df))){
    seed2<-as.numeric(subset_df[batch1,2])
    weight2<-as.numeric(subset_df[batch1,4:15])
    #samplesize2<-as.numeric(subset_df[batch1,1])
    dataframe2<-beta_arithematic_random_sequences_process(function1=sd,distype=dsnorm,parameters=c(location=0,scale=1),expect1=0,samplesize=samplesize2,seed1=seed2,weight1=weight2)
    dataframe_Process_Gaussian<-(cbind(dataframe2[2:13,2],weight2))
    all_Gaussian_sd<-rbind(all_Gaussian_sd,dataframe_Process_Gaussian)
  }
  
  Wmean00<-weighted.mean(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
  Wsd00<-weighted_SD(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
  Wmean02<-rbind(Wmean02,c(samplesize2,Wmean00))
  
  Wsd02<-rbind(Wsd02,c(samplesize2,Wsd00))
}

write.csv(Wmean02,file=paste("SD_beta-arith_random_Gaussian_240.csv", sep = ","), row.names=FALSE)
write.csv(Wsd02,file=paste("SD_beta-arith_random_Gaussian_240_sd.csv", sep = ","), row.names=FALSE)

simulatedbatch_Finitesample<-read.csv(file=paste("weights_norm3.csv", sep = ","))

Wmean03<-c()
Wsd03<-c()
for (samplesize2 in (5:100)){
  subset_df <- simulatedbatch_Finitesample[simulatedbatch_Finitesample$samplesize %in% samplesize2, ]
  all_Gaussian_sd<-c()
  dataframe2<-c()
  for(batch1 in (1:nrow(subset_df))){
    seed2<-as.numeric(subset_df[batch1,2])
    weight2<-as.numeric(subset_df[batch1,4:15])
    #samplesize2<-as.numeric(subset_df[batch1,1])
    dataframe2<-beta_arithematic_random_sequences_process(function1=sd,distype=dsnorm,parameters=c(location=0,scale=1),expect1=0,samplesize=samplesize2,seed1=seed2,weight1=weight2)
    dataframe_Process_Gaussian<-(cbind(dataframe2[2:13,2],weight2))
    all_Gaussian_sd<-rbind(all_Gaussian_sd,dataframe_Process_Gaussian)
  }
  
  Wmean00<-weighted.mean(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
  Wsd00<-weighted_SD(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
  Wmean03<-rbind(Wmean03,c(samplesize2,Wmean00))
  
  Wsd03<-rbind(Wsd03,c(samplesize2,Wsd00))
}
write.csv(Wmean03,file=paste("SD_beta-arith_random_Gaussian_360.csv", sep = ","), row.names=FALSE)
write.csv(Wsd03,file=paste("SD_beta-arith_random_Gaussian_360_sd.csv", sep = ","), row.names=FALSE)

simulatedbatch_Finitesample1<-read.csv(file=paste("weights_beta-arith_random.csv", sep = ","))

simulatedbatch_Finitesample2<-read.csv(file=paste("weights_norm1.csv", sep = ","))

simulatedbatch_Finitesample<-rbind(read.csv(file=paste("weights_norm3.csv", sep = ",")),as.data.frame(rbind(simulatedbatch_Finitesample1,simulatedbatch_Finitesample2)))


Wmean05<-c()
Wsd05<-c()
for (samplesize2 in (5:100)){
  subset_df <- simulatedbatch_Finitesample[simulatedbatch_Finitesample$samplesize %in% samplesize2, ]
  all_Gaussian_sd<-c()
  dataframe2<-c()
  for(batch1 in (1:nrow(subset_df))){
    seed2<-as.numeric(subset_df[batch1,2])
    weight2<-as.numeric(subset_df[batch1,4:15])
    #samplesize2<-as.numeric(subset_df[batch1,1])
    dataframe2<-beta_arithematic_random_sequences_process(function1=sd,distype=dsnorm,parameters=c(location=0,scale=1),expect1=0,samplesize=samplesize2,seed1=seed2,weight1=weight2)
    dataframe_Process_Gaussian<-(cbind(dataframe2[2:13,2],weight2))
    all_Gaussian_sd<-rbind(all_Gaussian_sd,dataframe_Process_Gaussian)
  }
  
  Wmean00<-weighted.mean(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
  Wsd00<-weighted_SD(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
  Wmean05<-rbind(Wmean05,c(samplesize2,Wmean00))
  
  Wsd05<-rbind(Wsd05,c(samplesize2,Wsd00))
}
write.csv(Wmean05,file=paste("SD_beta-arith_random_Gaussian_600.csv", sep = ","), row.names=FALSE)
write.csv(Wsd05,file=paste("SD_beta-arith_random_Gaussian_600_sd.csv", sep = ","), row.names=FALSE)

simulatedbatch_Finitesample1<-read.csv(file=paste("weights_beta-arith_random.csv", sep = ","))

simulatedbatch_Finitesample2<-read.csv(file=paste("weights_norm1.csv", sep = ","))

simulatedbatch_Finitesample_norm_5<-rbind(read.csv(file=paste("weights_norm3.csv", sep = ",")),as.data.frame(rbind(simulatedbatch_Finitesample1,simulatedbatch_Finitesample2)))

simulatedbatch_Finitesample_Weibull_1<-read.csv(file=paste("weights_Weibull_1.csv", sep = ","))

simulatedbatch_Finitesample_Weibull_2<-read.csv(file=paste("weights_Weibull_2.csv", sep = ","))

simulatedbatch_Finitesample_Weibull_5<-read.csv(file=paste("weights_Weibull_5.csv", sep = ","))

simulatedbatch_Finitesample_gnorm_2<-read.csv(file=paste("weights_gnorm2.csv", sep = ","))
simulatedbatch_Finitesample_gnorm_4<-read.csv(file=paste("weights_gnorm4.csv", sep = ","))

simulatedbatch_Finitesample_lnorm_025<-read.csv(file=paste("weights_lnorm025.csv", sep = ","))
simulatedbatch_Finitesample_lnorm_05<-read.csv(file=paste("weights_lnorm05.csv", sep = ","))

simulatedbatch_Finitesample_Pareto_7<-read.csv(file=paste("weights_Pareto_7.csv", sep = ","))
simulatedbatch_Finitesample_Pareto_10<-read.csv(file=paste("weights_Pareto_10.csv", sep = ","))
simulatedbatch_Finitesample_Pareto_15<-read.csv(file=paste("weights_Pareto_15.csv", sep = ","))


simulatedbatch_Finitesample_Weibull_gnorm<-as.data.frame(rbind(rbind(rbind(rbind(simulatedbatch_Finitesample_Weibull_1[1:288,1:15],simulatedbatch_Finitesample_Weibull_2[1:288,1:15]),simulatedbatch_Finitesample_Weibull_5[1:288,1:15]),simulatedbatch_Finitesample_gnorm_2[1:288,1:15]),simulatedbatch_Finitesample_gnorm_4[1:288,1:15]))

simulatedbatch_Finitesample_pareto_lnorm<-as.data.frame(rbind(rbind(rbind(rbind(simulatedbatch_Finitesample_Pareto_15[1:288,1:15],simulatedbatch_Finitesample_Pareto_10[1:288,1:15]),simulatedbatch_Finitesample_lnorm_025[1:288,1:15]),simulatedbatch_Finitesample_lnorm_05[1:288,1:15]),simulatedbatch_Finitesample_Pareto_7[1:288,1:15]))


simulatedbatch_Finitesample<-rbind(rbind(simulatedbatch_Finitesample_Weibull_gnorm[,1:15],simulatedbatch_Finitesample_pareto_lnorm[,1:15]),simulatedbatch_Finitesample2[1:1920,1:15])


Wmean5D<-c()
Wsd5D<-c()
for (samplesize2 in (5:100)){
  subset_df <- simulatedbatch_Finitesample[simulatedbatch_Finitesample$samplesize %in% samplesize2, ]
  all_Gaussian_sd<-c()
  dataframe2<-c()
  for(batch1 in (1:nrow(subset_df))){
    seed2<-as.numeric(subset_df[batch1,2])
    weight2<-as.numeric(subset_df[batch1,4:15])
    #samplesize2<-as.numeric(subset_df[batch1,1])
    dataframe2<-beta_arithematic_random_sequences_process(function1=sd,distype=dsnorm,parameters=c(location=0,scale=1),expect1=0,samplesize=samplesize2,seed1=seed2,weight1=weight2)
    dataframe_Process_Gaussian<-(cbind(dataframe2[2:13,2],weight2))
    all_Gaussian_sd<-rbind(all_Gaussian_sd,dataframe_Process_Gaussian)
  }
  
  Wmean00<-weighted.mean(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
  Wsd00<-weighted_SD(all_Gaussian_sd[,1],all_Gaussian_sd[,2])
  Wmean5D<-rbind(Wmean5D,c(samplesize2,Wmean00))
  
  Wsd5D<-rbind(Wsd5D,c(samplesize2,Wsd00))
}

Wmeanmonte_5<-c()
Wsdmonte1_5<-c()
for (samplesize2 in (5:100)){
  length2<-samplesize2
  sdnorm<-c()
  for (i in (1:600)){
    x<-rnorm(length2)
    sdnorm<-c(sdnorm,sd(x))
  }
  
  Wmeanmonte_5<-rbind(Wmeanmonte_5,c(samplesize2,mean(sdnorm)))
  Wsdmonte1_5<-rbind(Wsdmonte1_5,c(samplesize2,sd(sdnorm)))
}

write.csv(Wmeanmonte_5,file=paste("SD_Monte_Gaussian_600.csv", sep = ","), row.names=FALSE)
write.csv(Wsdmonte1_5,file=paste("SD_Monte_Gaussian_600_sd.csv", sep = ","), row.names=FALSE)

sqrt(mean((Wmeantrue[,2]-Wmeanmonte_5[,2])^2))
sqrt(mean((Wmeantrue[,2]-Wmean5D[,2])^2))
sqrt(mean((Wmeantrue[,2]-Wmeanarith[,2])^2))
sqrt(mean((Wmeantrue[,2]-Wmean1[,2])^2))
sqrt(mean((Wmeantrue[,2]-Wmeanmonte[,2])^2))
sqrt(mean((Wmeantrue[,2]-Wmean01[,2])^2))
sqrt(mean((Wmeantrue[,2]-Wmean02[,2])^2))
sqrt(mean((Wmeantrue[,2]-Wmean03[,2])^2))
sqrt(mean((Wmeantrue[,2]-Wmean05[,2])^2))



stopCluster(numCores)
registerDoSEQ()

