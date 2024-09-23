\name{coxphm}
\alias{coxphm}
\title{Time-to-Event Data Analysis with Missing Survival Times}
\description{Fits a pseudo Cox proprotional hazards model that allow us to analyze time-to-event data when survival times are missing for control groups.}
\usage{
coxphm(time, status, trt, z, beta0=NULL, time0=NULL, s=NULL, maxiter=1000, eps=0.01)
}
\arguments{
  \item{time}{Right-censored survival time (time is observed if trt=1. time is not observed if trt=0.)}
  \item{status}{Event indicator (status=1 if event, status=0 otherwise.)}
  \item{trt}{Treatment (or missing) indicator: trt=1 if treatment group (or no missing), trt=0 if control group (missing survival time).}
  \item{z}{Predictors (vector or matrix).}
  \item{beta0}{Initial value of regression parameters. (If beta0=NULL, estimated coefficient(s) from the logistic regression with status and z is used.)}
  \item{time0}{Initial value of (pseudo) survival times for trt=0. (If time=NULL, randomly selected time with replacement for trt=1 is used.)}
  \item{s}{Smoothed parameter. (If s=NULL, s=1/qnorm(1-n^(-2))*0.01 is used.)}
  \item{maxiter}{Number of maximim iteration.}
  \item{eps}{Stopping critiera.}
}
\details{Cox's proportional hazards model is not directly used to estimate a treatment effect when survival times for subjects in the control group(s) are missing. By regarding these missing survival times as nuisance parameters, the pseudo partial likelihood function is employed, which allows us to estimate the regression and nuisance parameters simultaneously with an unspecified baseline hazard function. In the pseudo partial likelihood, the smoothed parameter s is used to approximate risk sets as cumulative normal distributions with scale parameter s. Choosing a sufficient small s ensures that the pseudo partial likelihood is a good approximation of the partial likelihood. It is important to choose the initial value as close to the true value as possible. The estimated pseudo survival times range between 0 and Tmax, where Tmax is the maximum value of observed survival times for the treatment group.}

\value{
\item{conv}{Algorithm convergence: yes or not.}
\item{beta}{Estimated regression parameter(s): beta: estimated coefficient, se: standard error; lcl: lowr confidence limit, ucl: upper confidence limit, statistics: test statistics; pvalue: pvalue.}
\item{eta}{Estimated pseudo survival time.}
\item{loglik}{Log pseudo-partial-likelihood value.}
\item{iter}{Number of iterations.}
}

\references{Proportional hazards model when time-origin is not identifiable for control group (in-progress)}
\author{Yunro Chung [aut, cre]}

\examples{
#Mayo's pbc dataset from the survival package.
pbc1=pbc[1:200,] #first 200 patients

time=pbc1$time
status=pbc1$status
status[which(status==1)]=0 #transplant
status[which(status==2)]=1 #death

trt=pbc1$trt
trt[which(trt==2)]=0 #0 for placebo, 1 for treatment

age=pbc1$age
z=cbind(trt,age)
colnames(z)=c("trt","age")

#0. Cox model
fit0=coxph(Surv(time,status)~z)

#1. Pseudo-Cox model
#1.1. initial value
beta0=fit0$coefficients
time0=time[which(trt==0)]

#1.2. Survival times are missing if trt=0
time[which(trt==0)]=NA

#1.3. fits pseudo-Cox
fit1=coxphm(time, status, trt=trt, z=z, beta0=beta0, time0=time0)

#2. Cox vs pseudo-cox (almost identifical)
print(summary(fit0)$coefficient)
print(fit1$beta)
print(time0-fit1$eta)

#3. Subsequent analyses after fitting pseudo-cox
time[which(trt==0)]=fit1$eta
survfit(Surv(time,status)~trt) #Kaplan-Meiere
survdiff(Surv(time,status)~trt)#Log-rank test
}

