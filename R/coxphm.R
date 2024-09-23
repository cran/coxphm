coxphm=function(time, status, trt=NULL, z, beta0=NULL, time0=NULL, s=NULL, maxiter=1000, eps=0.01){

  #0. set up
  time.min=0 #min(time,na.rm=TRUE)
  time.max=max(time,na.rm=TRUE)

  n=length(time)
  if(is.null(s))
    s=1/qnorm(1-n^(-2))*0.001
  if(is.null(beta0)){
    zmat=as.matrix(z,nrow=n)
    beta0=coef(glm(status ~ zmat, family = "binomial"))[-1] #from logistic
  }

  if(is.null(trt)){
    trt=rep(NA,n)
    trt[is.na(time)]=0
    trt[!is.na(time)]=1
  }

  if(is.null(time0)){
    time0=sample(time[which(trt==1)],length(time[which(trt==0)]),replace=TRUE) #random selection
  }

  #1. define subgroups and variables
  x=as.numeric(!is.na(time)) #missing indicator should be same as x=1 for treatment; x=0 for control (same as trt)
  df=data.frame(time=time,status=status,x=x,z=z)
  dfC=df[which(x==0),]
  dfT=df[which(x==1),]

  n=length(time)
  mC=nrow(dfC)
  mT=nrow(dfT)
  p=ncol(dfC)-3 #time, status, x

  if(p==1){
    namez=names(z)
  }else{
    namez=colnames(z)
  }

  #2. define new data structure
  dfC0=list(m=mC,time=dfC$time, #NA for time
            status=dfC$status,wstatus=which(dfC$status==1),
            x=dfC$x,
            z=as.matrix(dfC[,-c(1,2,3)],nrow=mC),  #-1,2,3 for time,status,x
            zb=NA,ezb=NA
  )

  dfT0=list(m=mT,time=dfT$time,
            status=dfT$status,wstatus=which(dfT$status==1),
            x=dfT$x,
            z=as.matrix(dfT[,-c(1,2,3)],nrow=mT),
            zb=NA,ezb=NA)

  #3.initial value
  beta0=matrix(beta0,ncol=1)

  dfC0$time=time0

  dfC0$zb=c(dfC0$z %*% beta0)
  dfC0$ezb=exp(dfC0$zb)

  dfT0$zb=c(dfT0$z %*% beta0)
  dfT0$ezb=exp(dfT0$zb)

  #3. NR + optim
  conv="no"
  dist=1
  iter=0
  eta0=time0

  while(conv=="no"){
    iter=iter+1

    #1. update beta given eta using newton-raphson
    drvt1.res=drvt1(dfC=dfC0, dfT=dfT0, p=p, s=s)
    try1=try(beta1<-beta0-ginv(drvt1.res$H)%*%(drvt1.res$U),silent=TRUE)
    if(class(try1)[1]=="try-error")
      break

    dfC1=dfC0
    dfC1$zb=c(dfC1$z %*% beta1)
    dfC1$ezb=exp(dfC1$zb)

    dfT1=dfT0
    dfT1$zb=c(dfT1$z %*% beta1)
    dfT1$ezb=exp(dfT1$zb)

    #2. update eta given beta using optim
    par=c(eta0)
    TT=data.frame(time=dfC1$time,status=dfC1$status,z=dfC1$z, trt=dfC1$x, s=s, zb=dfC1$zb, ezb=dfC1$ezb)
    CC=data.frame(time=dfT1$time,status=dfT1$status,z=dfT1$z, trt=dfT1$x, s=s, zb=dfT1$zb, ezb=dfT1$ezb)
    data=rbind(TT,CC)

    #optim.fit <- optim(par = par, fn = n_ps_pll, data=data, control=list(maxit=maxiter), method = "BFGS", hessian = FALSE)
    try2=try(optim.fit <- optim(par = par, fn = n_ps_pll, data=data, control=list(maxit=maxiter), method = "L-BFGS-B", hessian = FALSE, lower=time.min, upper=time.max),silent=TRUE)
    if(class(try2)[1]=="try-error")
      break

    eta1=optim.fit$par

    #3. euclidean distance
    dist=sqrt(sum((beta1-beta0)^2))+sqrt(sum((eta1-eta0)^2))


    if(is.na(dist))
      break

    if(dist<eps)
      conv="yes"

    dfC0=dfC1
    dfT0=dfT1

    beta0=beta1
    eta0=eta1

    #print(round(c(iter,dist),3))
    if(iter>maxiter)
      break
  }

  if(conv=="no"){ #not converged
    res=list(conv="no",beta=NA,eta=NA,loglik=NA)
  }else{
    qz=qnorm(1 - 0.05/2)

    beta=beta1
    eta=eta1

    drvt1.res=drvt1(dfC=dfC1, dfT=dfT1, p=p, s=s)
    se=sqrt(-diag(ginv(drvt1.res$H)))
    lcl=beta-qz*se
    ucl=beta+qz*se

    beta.mat=data.frame(beta=beta,se=se,lcl=lcl,ucl=ucl,statistics=beta/se)
    beta.mat$pvalue=(1-pnorm(abs(beta.mat$statistics)))*2

    rownames(beta.mat)=namez

    loglik=lpl(dfC1,dfT1,s)

    res=list(conv="yes",beta=beta.mat,eta=eta,loglik=loglik,iter=iter)
  }

  return(res)
}
