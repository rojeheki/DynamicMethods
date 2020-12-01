# Based on "Estimating a nonlinear rational expectations commodity price model
# with unobservable state variables" by Angus Deaton and Guy Laroque
# Journal of Applied Econometrics, Vol. 10, S9-S40 (1995)
# https://doi.org/10.1002/jae.3950100503

library(cubature)
library(ggplot2)

# 3.1 Discretization of the Harvest

generateZTransitionMatrix = function() {
  f = function(x) {
    exp(-(x[1]^2-2*s$rho*x[1]*x[2]+x[2]^2)/(2*(1-s$rho^2)))/(2*pi*sqrt(1-s$rho^2))
  }
  theta = qnorm(0:s$n/s$n)
  s$Tz = matrix(nrow = s$n, ncol = s$n)
  for(i in 1:s$n) {
    for(j in 1:i){
      s$Tz[i,j]=s$Tz[j,i]=adaptIntegrate(f, lowerLimit=c(theta[i],theta[j]),upperLimit=c(theta[i+1],theta[j+1]))$integral*s$n
    }
  }
}

calculateZLevels = function(){
  theta = qnorm(0:s$n/s$n)
  s$Z = (-(dnorm(theta[-1])-dnorm(theta[-(s$n+1)]))*s$n)
}

s = new.env()

setZDiscretization = function(n,rho){
  s$n=n
  s$rho=rho
  generateZTransitionMatrix()
  calculateZLevels()
}

setZDiscretization(10,0.7)

expectedHarvest = function(currentZLvl){
  return (sum(s$Tz[currentZLvl,]*s$Z))
}

varianceExpectedHarvest = function(currentZLvl){
  mean = expectedHarvest(currentZLvl)
  return (sum(s$Tz[currentZLvl,]*(s$Z)^2)-mean^2)
}

# Figure 1
setZDiscretization(10,0.7)
e = c()
for (i in 1:10){
  e[i]=expectedHarvest(i)
}

plot(-2:2, -2:2*0.7, type="l")
points(s$Z[1:10],e, type="l")

# Figure 2 - Adjustment of the actual conditional variance by 1/(1-rho^2) necessary to replicate the chart
plot(c(-2,2),c(0,1.5))
for (rho in c(0.1,0.3,0.5,0.7,0.9)){
  setZDiscretization(10,rho)
  v = c()
  for (i in 1:10){
    v[i]=varianceExpectedHarvest(i)/(1-rho^2)
  }
  points(s$Z[1:10],v, type="l")
}

# Figure 3
setZDiscretization(10,0.7)
lvls = c(5)
r = runif(200)
for (i in 1:200) {
  j=1
  while (r[i] > s$Tz[lvls[i],j]) {
    r[i] = r[i]-s$Tz[lvls[i],j]
    j = j+1
  }
  lvls[i+1] = j
}
plot(1:201, s$Z[lvls], type="l")

# (Inverse) Demand function

s$P = function(x){return (s$a+s$b*x)}
s$D = function(p){return ((p-s$a)/s$b)}

# set static parameters of the QMLE

q = new.env()
q$constants["r"] = 0.05

getR = function(){
  return (q$constants["r"])
}

# Transform to and from generic parameter vector theta using formula (43) (not the (43) that's actually (42))

toTheta = function(a,b,delta){
  s$a = a
  s$b = b
  s$delta = delta
  q$theta = c(a, log(-b), log(delta+0.05))
}

fromTheta = function(theta){
  q$theta = theta
  s$a = theta[1]
  s$b = -exp(theta[2])
  s$delta = -0.05+exp(theta[3])
}

toTheta(0.2, -0.15, 0.12)


# 3.2 Discretization of the availability

setXDiscretization = function (m) {
  # minimum is 0 stored + minimum harvest, maximum is maximum harvest/delta with a safety factor allowing for delta to shrink
  safety_factor = 2
  s$X = (0:(m-1))/m*(s$Z[s$n]/s$delta*safety_factor-s$Z[1])+s$Z[1]
  s$m=m
}

setXDiscretization(20)

# 3.2 Iterate to improve price function f
# f(x,z) is stored as a matrix with f[k,i] representing f(X[k],Z[i])

initF = function(){
  s$f = matrix(nrow=s$m, ncol=s$n)
  for (k in 1:s$m){
    s$f[k,]=max(s$P(s$X[k]), 0)
  }
}


# Using formula (28)

iterateFOnce = function (){
  newF = matrix(nrow=s$m, ncol=s$n)
  beta = (1-s$delta)/(1+getR())
  splines=list()
  for (i in 1:s$n){
    splines=append(splines,splinefun(s$X, s$f[,i]))
  }
  futureValue = function(x,i){
    v = 0
    y = s$D(splines[[i]](x))
    for (j in 1:s$n){
      v = v + s$Tz[i,j]*splines[[j]](s$Z[j]+(1-s$delta)*(x-y))
    }
    return (beta*v)
  }
  for (k in 1:s$m){
    currentPrice=s$P(s$X[k])
    for (i in 1:s$n){
      newF[k,i]=max(futureValue(s$X[k],i),currentPrice)
    }
  }
  s$f=newF
}

iterateF = function(dCutoff = 1e-5, nCutoff = 100) {
  initF()
  lastF = s$f+Inf
  for (i in 1:nCutoff) {
    deviation = max(abs(s$f - lastF))
    if(deviation<dCutoff) {
      print(i)
      return (TRUE)
    }
    lastF = s$f
    iterateFOnce()
  }
  return (FALSE)
}


iterateF()


plot(s$X,s$f[,1],type="l")
for(i in 2:s$n){
  points(s$X,s$f[,i],type="l")
}



# 3.3 Generate the combined storage/harvest transition matrix

# calculates g(x,z) from formulas (36), (39)
expectedStorage = function(currentXLvl, currentZLvl) {
  return ((1-s$delta)*(s$X[currentXLvl]-s$D(s$f[currentXLvl,currentZLvl]))+s$rho*s$Z[currentZLvl])
}


# XZ Transition matrix will be ordered by Z on a large scale, X on a small scale
# X1Z1, X2Z1, ..., XmZ1, X1Z2, X2Z2, ..., xmZn

generateXZTransitionMatrix = function(){
  theta = qnorm(0:s$n/s$n)
  s$Txz = matrix(nrow = s$m*s$n, ncol = s$m*s$n)
  xd = append(append(Inf,(s$X[2:s$m]-s$X[1:(s$m-1)])/2),Inf) # delta/2 between X, xd[i] is left of X[i], xd[i+1] right
  for (i in 1:s$m) {
    for (j in 1:s$n) {
      c_row = (i-1)*s$n + j
      for (k in 1:s$m) {
        for (l in 1:s$n) {
          c_column = (k-1)*s$n + l
          xdc = s$X[i] - expectedStorage(k,l)
          critXL = xdc - xd[i]
          critXR = xdc + xd[i+1]
          critZL = theta[j]-s$rho*s$Z[l]
          critZR = theta[j+1]-s$rho*s$Z[l]
          critL = max(critXL,critZL)
          critR = min(critXR,critZR)
          if (critL > critR){
            s$Txz[c_row,c_column]=0
          } else {
            s$Txz[c_row,c_column]=pnorm(critR)-pnorm(critL)
          }
        }
      }
    }
  }
}

setZDiscretization(10,0.7)
setXDiscretization(20)
iterateF()


generateXZTransitionMatrix()


# Calculate invariant distribution

s$id = eigen(s$Txz)$vectors[,1]

# Visualize invariant distribution

idVis = data.frame()
for (i in 1:s$m){
  for (j in 1:s$n){
    nl = data.frame(x=s$X[i], z=s$Z[j], p=as.numeric(s$id[(i-1)*s$n+j]))
    idVis = rbind(idVis,nl)
  }
}


idPlot = ggplot(idVis, aes(x=z,y=x,z=p)) + geom_contour_filled()
idPlot


# Generate Monte-Carlo data

s$t = 100
s$mc = data.frame(z = 1, x = 1, p=s$f[5,10])
r = runif(s$t)
for (i in 1:(s$t-1)) {
  j=1
  while (r[i] > s$Txz[j,(s$mc[i,1]-1)*s$n+s$mc[i,2]]) {
    r[i] = r[i]-s$Txz[j,(s$mc[i,1]-1)*s$n+s$mc[i,2]]
    j = j+1
  }
  s$mc[i+1,1] = (j-1)%%s$n+1
  s$mc[i+1,2] = floor((j-1)/s$n)+1
  s$mc[i+1,3] = s$f[s$mc[i+1,1],s$mc[i+1,2]]
}

plot(1:s$t,s$mc$p, type="l")

# 5.1 Estimate (inverse) price function for AR case
# Requires data



# 5.1 Calculate 1-period-ahead expectations and variances using formulas (45) and (46)

# Can the iid case simply use the AR functions for expectation and variance or does it  
# have to be calculated using a separate procedure?

