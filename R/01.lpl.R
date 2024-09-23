#log partial likelihood
lpl=function(dfC, dfT,s){
  loglik=sum(dfC$zb[dfC$wstatus])+sum(dfT$zb[dfT$wstatus])
  for(i in dfC$wstatus)
    loglik=loglik-log(sum(I1(dfC$time[i],dfC$time,s)*dfC$ezb)+sum(I1(dfC$time[i],dfT$time,s)*dfT$ezb))
  for(i in dfT$wstatus)
    loglik=loglik-log(sum(I1(dfT$time[i],dfC$time,s)*dfC$ezb)+sum(I0(dfT$time[i],dfT$time)*dfT$ezb))
  return(loglik)
}
