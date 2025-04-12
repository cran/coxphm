\name{coxphm}
\alias{coxphm}
\title{Time-to-Event Data Analysis with Missing Survival Times}
\description{Fits a pseudo Cox proprotional hazards model that allow us to analyze time-to-event data when survival times are missing for control groups.}
\usage{
coxphm(time, status, trt, z, beta0, time0, Atime, Btime, u, s, maxiter, eps)
}
\arguments{
  \item{time}{Right-censored survival time (time is observed if trt=1. time is not observed if trt=0.)}
  \item{status}{Event indicator (status=1 if event, status=0 otherwise.)}
  \item{trt}{Treatment (or missing) indicator: trt=1 if treatment group (or no missing), trt=0 if control group (missing survival time).}
  \item{z}{Predictors (vector or matrix), where z[,1] must be the same as trt.}
  \item{beta0}{Initial value of regression parameters.}
  \item{time0}{Initial value of (pseudo) survival times for trt=0.}
  \item{Atime}{Duration from treatment available to time-origin.}
  \item{Btime}{Duration from treatment available to survival time. Missing if trt=1.}
  \item{u}{A variable to estimate time0. See the below.}
  \item{s}{Smoothed parameter. (If s=NULL, s=1/qnorm(exp(-n^(-2)))*0.01 is used, where n is the number of samples.)}
  \item{maxiter}{Number of maximim iteration. Default: maxiter=1000.}
  \item{eps}{Stopping critiera. Default: eps=0.01.}
}
\details{Cox's proportional hazards model is not directly used to estimate treatment effects when survival times for subjects in the control group are missing. By regarding these missing survival times as nuisance parameters, the pseudo partial likelihood function is used to estimate the regression and nuisance parameters simultaneously with an unspecified baseline hazard function.

In the pseudo partial likelihood, the smoothed parameter s is used to approximate risk sets as cumulative normal distributions with s as the scale parameter. Choosing a sufficiently small s ensures that the pseudo partial likelihood accurately approximates the partial likelihood.

The method depends on the choice of initial value, and hence, it is crucial to choose the initial value as close to the true value as possible. If this is not possible, a data-driven initial value can is utilized as follows. If beta0=NULL, a logistic regression is used to model the probability of status=1 using z as a predictor, and beta0 is set to the regression coefficient from the logistic regression. If time0=NULL, a linear regression is used to model Atime using u as a predictor, and time0 is set to Btime minus estimated Atime. Here, u and Btime are observed for all samples, whereas Atime is observed for trt=1. Here, Atime,Btime and u are not used if time0 is input.}

\value{
\item{conv}{Algorithm convergence: yes or not.}
\item{beta}{Estimated regression parameter: beta: estimated coefficient, se: standard error; lcl: 95\% lower confidence limit, ucl: 95\% upper confidence limit, statistics: test statistics; pvalue: pvalue.}
\item{eta}{Estimated pseudo survival time.}
\item{loglik}{Log pseudo-partial-likelihood value.}
\item{iter}{Number of iterations.}
}

\references{Pseudo partial likelihood method for proportional hazards models when time origin is missing for control group with applications to SARS-CoV-2 Seroprevalence Study (in-progress)}
\author{Yunro Chung [aut, cre]}

\examples{
#Mayo's pbc dataset from the survival package.
pbc1=pbc[1:200,] #use the first 200 patients

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
#1.1. set (true) initial value
beta0=fit0$coefficients
time0=time[which(trt==0)]

#1.2. Survival times are missing if trt=0
time[which(trt==0)]=NA

#1.3. fits pseudo-Cox
fit1=coxphm(time, status, trt=trt, z=z, beta0=beta0, time0=time0)

#2. Cox vs pseudo-cox
print(summary(fit0)$coefficient)
print(fit1$beta)
}
