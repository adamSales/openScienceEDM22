library(tidyverse)
ns=c(50,100,500,1000)

### covariate selection

makeDat <- function(n,p){

  ## 10 covariates, all normally distributed
  dat <- as.data.frame(scale(matrix(rnorm(n*p),ncol=p)))
  dat$treatment <- rep(c(0,1),n/2)
  dat$Y <- rnorm(n) ## no trt effect. R2=.5

  dat
}

backwards <- function(dat){
  mod0 <- lm(Y~.^2,data=dat)
  modFinal<-
    step(mod0,
         scope=list(upper=Y~.^2,lower=Y~treatment),
         direction = 'backward',
         trace=0)
  summary(modFinal)$coef['treatment',c('Estimate','Pr(>|t|)')]
}

oracle <- function(dat){
  mod <- lm(Y~treatment+V1,data=dat)
  summary(mod)$coef['treatment',c('Estimate','Pr(>|t|)')]
}

ttest <- function(dat){
  mod <- lm(Y~treatment,data=dat)
  summary(mod)$coef['treatment',c('Estimate','Pr(>|t|)')]
}

doAll <- function(i=1,n,p){
  dat <- makeDat(n=n,p=p)
  rbind(backwards=backwards(dat),
#    oracle=oracle(dat),
    ttest=ttest(dat))
}



library(parallel)
cl <- makeCluster(12)
clusterExport(cl,list('makeDat','backwards','oracle','ttest','doAll'))
design <- expand.grid(n=ns,p=c(2,5,10))
design <- subset(design,n> (p+1)*(p+2)/2+1)
res <- lapply(1:nrow(design),
              function(i){
                  cat(round(i/nrow(design)*100),'% ',sep='')
                  parLapply(cl,1:1000,doAll,n=design$n[i],p=design$p[i])
              })
cat('\n')
save(res,design,file='variableSelection.RData')

lev <- sapply(res, function(run) rowMeans(sapply(run,function(x) x[,2]<0.05)))
lev=cbind(design,t(lev))
save(res,design,lev,file='variableSelection.RData')

rej=sapply(res, function(run) rowSums(sapply(run,function(x) x[,2]<0.05)))

backTest=lapply(lev$backwards*1000,binom.test,n=1000)
ttestTest=lapply(lev$ttest*1000,binom.test,n=1000)

backCI <- sapply(backTest,function(x) x$conf.int)
ttestCI <- sapply(backTest,function(x) x$conf.int)

lev=cbind(lev,backCI=t(backCI))

pdat=lev%>%
    pivot_longer(c(backwards,ttest),names_to=c('method'),values_to='error')

ci=sapply(pdat$error,function(x) binom.test(x*1000,1000)$conf.int)
pdat=cbind(pdat,ci=t(ci))

pdat$n=pdat$n+ifelse(pdat$method=='ttest',10,-10)

pdat%>%
    ggplot(aes(n,error,ymin=ci.1,ymax=ci.2,color=method,group=method))+geom_point()+geom_line()+
    geom_errorbar(width=0)+
    geom_hline(yintercept=0.05,linetype='dotted')+
    facet_wrap(~p)
ggsave('covSelection.jpg')



################ log transformation
skewness <- function(x){
  n <- length(x)
  mu <- mean(x)

  (n^2/((n-1)*(n-2))*mean((x-mu)^3)/(sd(x)^(3/2)))
}

makeDat <- function(n)
  data.frame(Y=exp(rnorm(n)),treatment=rep(c(0,1),n/2))

analyze <- function(dat){
  skew <- max(skewness(dat$Y[dat$treatment==1]),skewness(dat$Y[dat$treatment==0]))
  ttest <- summary(lm(Y~treatment,data=dat))$coef['treatment',]
  ttestLog <- summary(lm(log(Y)~treatment,data=dat))$coef['treatment',]
  chooseSkew <- ifelse(skew>6,ttestLog['Pr(>|t|)'],ttest['Pr(>|t|)'])
  chooseP <- min(ttestLog['Pr(>|t|)'],ttest['Pr(>|t|)'])
  c(lin=ttest['Pr(>|t|)'],log=ttestLog['Pr(>|t|)'],choose=chooseSkew,chooseP=chooseP,skew)
}

clusterExport(cl,list('makeDat','analyze','skewness'))

res <- lapply(ns,
              function(n) parLapply(cl,1:10000,function(i) analyze(makeDat(n))))
save(res,ns,file='logTransformation.RData')
lev=sapply(res, function(x) colMeans(do.call('rbind',x)[,1:4]<0.05))
save(res,ns,lev,file='logTransformation.RData')

###################### excluding outliers

outlier <- function(x){
  quants <- quantile(x,c(0.25,0.75))
  iqr <- quants[2]-quants[1]
  (x< quants[1]-1.5*iqr)|(x> quants[2]+1.5*iqr)
}

sim <- function(n){
  Y <- rt(n,3)
  treatment <- rep(c(0,1),each=n/2)
  Yt <- Y[1:(n/2)]
  Yc <- Y[(n/2+1):n]
  Ytd <- Yt[!outlier(Yt)]
  Ycd <- Yc[!outlier(Yc)]

  c(
    ttest=t.test(Yt,Yc)$p.value,
    dropOut=t.test(Ytd,Ycd)$p.value)
}

clusterExport(cl,list('outlier','sim'))

res <- lapply(ns,
              function(n) parLapply(cl,1:10000,function(i) sim(n=n)))
save(res,ns,file='excludeOutliers.RData')

lev=sapply(res, function(x) colMeans(do.call('rbind',x)<0.05))
save(res,ns,lev,file='logTransformation.RData')


stopCluster(cl)
