#1. negative partial log-likelihood function for optim algo
n_ps_pll = function(par,data){ #s is scaling parameter for ind.approx

  #2.1. Initial param values needed for optimization
  eta  = par

  #2.2 Replace T0 for the unvaccinated group with the eta parms
  s=data$s[1] #smoothed parameter
  data$y=data$time
  data$y[which(data$trt==0)]=eta

  #2.3. likelihood
  dfT=data[which(data$trt==1),] #trt, vaccinated
  dfC=data[which(data$trt==0),] #ctl, unvaccinated

  wstatusT=which(dfT$status==1) #indicator at event time
  wstatusC=which(dfC$status==1)

  loglik=sum(data$status*data$zb)
  for(i in wstatusC)
    loglik=loglik-log(sum(I1(dfC$y[i],dfC$y,s)*dfC$ezb)+sum(I1(dfC$y[i],dfT$y,s)*dfT$ezb))
  for(i in wstatusT)
    loglik=loglik-log(sum(I1(dfT$y[i],dfC$y,s)*dfC$ezb)+sum(I0(dfT$y[i],dfT$y)*dfT$ezb))

  nloglik=-(loglik)
  return(nloglik)
}
