#drvt for beta
drvt1=function(dfC, dfT, p, s){
  if(p==1){; U=sum(dfC$z[dfC$wstatus,])+sum(dfT$z[dfT$wstatus,])
  }else{;    
    if(length(dfC$wstatus)==1){
      U1=dfC$z[dfC$wstatus,]
    }else{
      U1=apply(dfC$z[dfC$wstatus,],2,sum)
    }
    if(length(dfT$wstatus)==1){
      U2=dfT$z[dfT$wstatus,]
    }else{
      U2=apply(dfT$z[dfT$wstatus,],2,sum)
    }
    U=matrix(U1+U2,nrow=p)
  }
  H=matrix(0,nrow=p,ncol=p)
  for(i in dfC$wstatus){
    den.i=sum(I1(dfC$time[i],dfC$time,s)*dfC$ezb)+sum(I1(dfC$time[i],dfT$time,s)*dfT$ezb)
    num1.i=0
    num2.i=0
    for(j in 1:dfC$m){
      num1.i=num1.i+I1(dfC$time[i],dfC$time[j],s)*dfC$ezb[j]*dfC$z[j,]
      num2.i=num2.i+I1(dfC$time[i],dfC$time[j],s)*dfC$ezb[j]*(dfC$z[j,]%*%t(dfC$z[j,]))
    }
    for(j in 1:dfT$m){
      num1.i=num1.i+I1(dfC$time[i],dfT$time[j],s)*dfT$ezb[j]*dfT$z[j,]
      num2.i=num2.i+I1(dfC$time[i],dfT$time[j],s)*dfT$ezb[j]*(dfT$z[j,])%*%t(dfT$z[j,])
    }
    U=U-num1.i/den.i
    H=H-num2.i/den.i+(num1.i%*%t(num1.i))/(den.i^2)
  }
  for(i in dfT$wstatus){
    den.i=sum(I1(dfT$time[i],dfC$time,s)*dfC$ezb)+sum(I0(dfT$time[i],dfT$time)*dfT$ezb)
    num1.i=matrix(0,nrow=p)
    num2.i=matrix(0,nrow=p,ncol=p)
    for(j in 1:dfC$m){
      num1.i=num1.i+I1(dfT$time[i],dfC$time[j],s)*dfC$ezb[j]*dfC$z[j,]
      num2.i=num2.i+I1(dfT$time[i],dfC$time[j],s)*dfC$ezb[j]*(dfC$z[j,])%*%t(dfC$z[j,])
    }
    for(j in 1:dfT$m){
      num1.i=num1.i+I0(dfT$time[i],dfT$time[j])*dfT$ezb[j]*dfT$z[j,]
      num2.i=num2.i+I0(dfT$time[i],dfT$time[j])*dfT$ezb[j]*(dfT$z[j,])%*%t(dfT$z[j,])
    }
    U=U-num1.i/den.i
    H=H-num2.i/den.i+(num1.i%*%t(num1.i))/(den.i^2)
  }
  return(list(U=U,H=H))
}
