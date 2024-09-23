Phi=function(x,s)
  pnorm(x/s)

I0=function(i,j)   I(j>=i) #equivalent to I(j-i>=0)
I1=function(i,j,s) Phi(j-i,s)

phi=function(x,s)
  dnorm(x/s)
I2=function(i,j,s) phi(j-i,s)/s

#I1=function(i,j,s) Phi(j-i,s)*I(i!=j)  #I(i!=j) is used because I(i>=j) is 0, and its smooth approximation is zero.
#I2=function(i,j,s) phi(j-i,s)/s*I(i!=j)

#  dnorm(x/s)
#I3=function(i,j,s) phi(j-i,s)*(j-i)/(s^3) #for the second derivative of eta
