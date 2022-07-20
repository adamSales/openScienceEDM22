### covariate selection

makeDat <- function(n){

  ## 10 covariates, all normally distributed
  dat <- as.data.frame(scale(matrix(rnorm(n*10),ncol=10)))
  dat$treatment <- rep(c(0,1),n/2)
  dat$Y <- rnorm(n) ## no trt effect. R2=.5

  dat
}

backwards <- function(dat){
  mod0=lm(Y~.^2,data=dat)
  modFinal<-
    step(mod0,
         scope=list(upper=Y~.^2,lower=Y~treatment),
         direction = 'backward',
         trace=0)
  summary(modFinal)$coef['treatment',c('Estimate','Pr(>|t|)')]
}

oracle <- function(dat){
  mod=lm(Y~treatment+V1,data=dat)
  summary(mod)$coef['treatment',c('Estimate','Pr(>|t|)')]
}

ttest <- function(dat){
  mod=lm(Y~treatment,data=dat)
  summary(mod)$coef['treatment',c('Estimate','Pr(>|t|)')]
}

doAll=function(i=1,n){
  dat <- makeDat(n=n)
  rbind(backwards=backwards(dat),
#    oracle=oracle(dat),
    ttest=ttest(dat))
}

library(parallel)
cl=makeCluster(4)
clusterExport(cl,list('makeDat','backwards','oracle','ttest','doAll'))
system.time(res <- parLapply(cl,1:1000,doAll,n=200))
lev=sapply(res,function(x) x[,2]<0.05)
rowMeans(lev)


################
skewness <- function(x){
  n=length(x)
  mu=mean(x)

  (n^2/((n-1)*(n-2))*mean((x-mu)^3)/(sd(x)^(3/2)))
}

makeDat=function(n)
  data.frame(Y=exp(rnorm(n)),treatment=rep(c(0,1),n/2))

analyze=function(dat){
  skew=skewness(dat$
  ttest=t.test(Y~treatment,data=dat)$p.value
  ttestLog=t.test(log(Y)~treatment,data=dat)$p.value
  choose=ifelse(max(skewness(dat$Y[dat$treatment==1]),skewness(dat$Y[dat$treatment==0]))>6,ttestLog,ttest)
  c(ttest,ttestLog,choose)
}

system.time(res<-replicate(10000,analyze(makeDat(50))))


###################### excluding outliers

outlier=function(x){
  iqr=IQR(x)
  quants=quantile(x,c(0.25,0.75))
  iqr=quants[2]-quants[1]
  (x< quants[1]-1.5*iqr)|(x> quants[2]+1.5*iqr)
}

sim=function(n){
  Y=rt(n,3)
  treatment=rep(c(0,1),each=n/2)
  Yt=Y[1:(n/2)]
  Yc=Y[(n/2+1):n]
  Ytd=Yt[!outlier(Yt)]
  Ycd=Yc[!outlier(Yc)]

  c(
    ttest=t.test(Yt,Yc)$p.value,
    dropOut=t.test(Ytd,Ycd)$p.value)
}
