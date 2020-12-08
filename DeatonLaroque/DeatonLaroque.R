# Based on Angus Deaton and Guy Laroque. Estimating a nonlinear rational expectations commodity price model with unobservable state variables
# Journal of Applied Econometrics, Vol. 10, S9-S40 (1995)
# https://doi.org/10.1002/jae.3950100503

library(cubature)
library(ggplot2)

# Import and organize data

dat = read.csv2("TimeSeries.csv", sep=",")
com = list(cocoa=new.env(), coffee=new.env(), copper=new.env(), rice=new.env(), sugar=new.env(), tin=new.env())

com$cocoa$pDat = c(dat$P_Cocoa)/c(dat$US_CPI)
com$cocoa$years = c(dat$Year)

com$coffee$pDat = c(dat$P_Coffee)/c(dat$US_CPI)
com$coffee$years = c(dat$Year)

com$copper$pDat = c(dat$P_Copper)/c(dat$US_CPI)
com$copper$years = c(dat$Year)

com$rice$pDat = c(dat$P_Rice)/c(dat$US_CPI)
com$rice$years = c(dat$Year)

com$sugar$pDat = c(dat$P_Sugar)/c(dat$US_CPI)
com$sugar$years = c(dat$Year)

com$tin$pDat = c(dat$P_Tin)/c(dat$US_CPI)
com$tin$years = c(dat$Year)

calculateT = function(e){
  e$t = length(e$years)
}

for (e in com) {
  calculateT(e)
}


# 3.1 Discretization of the Harvest

generateZTransitionMatrix = function(e) {
  f = function(x) {
    exp(-(x[1]^2-2*e$rho*x[1]*x[2]+x[2]^2)/(2*(1-e$rho^2)))/(2*pi*sqrt(1-e$rho^2))
  }
  theta = qnorm(0:e$n/e$n)
  e$Tz = matrix(nrow = e$n, ncol = e$n)
  for(i in 1:e$n) {
    for(j in 1:i){
      e$Tz[i,j]=e$Tz[j,i]=adaptIntegrate(f, lowerLimit=c(theta[i],theta[j]),upperLimit=c(theta[i+1],theta[j+1]))$integral*e$n
    }
  }
}

calculateZLevels = function(e){
  theta = qnorm(0:e$n/e$n)
  e$Z = (-(dnorm(theta[-1])-dnorm(theta[-(e$n+1)]))*e$n)
}

setZDiscretization = function(n,rho,e){
  e$n=n
  e$rho=rho
  generateZTransitionMatrix(e)
  calculateZLevels(e)
}

setZDiscretization(10,0.7,com$cocoa)

expectedHarvest = function(currentZLvl, e){
  return (sum(e$Tz[currentZLvl,]*e$Z))
}

varianceExpectedHarvest = function(currentZLvl, e){
  mean = expectedHarvest(currentZLvl, e)
  return (sum(e$Tz[currentZLvl,]*(e$Z)^2)-mean^2)
}

# Figure 1
setZDiscretization(10,0.7, com$cocoa)
eH = c()
for (i in 1:10){
  eH[i]=expectedHarvest(i, com$cocoa)
}

plot(-2:2, -2:2*0.7, type="l")
points(com$cocoa$Z[1:10],eH, type="l")

# Figure 2 - Adjustment of the actual conditional variance by 1/(1-rho^2) necessary to replicate the chart
plot(c(-2,2),c(0,1.5))
for (rho in c(0.1,0.3,0.5,0.7,0.9)){
  setZDiscretization(10,rho,com$cocoa)
  v = c()
  for (i in 1:10){
    v[i]=varianceExpectedHarvest(i,com$cocoa)/(1-rho^2)
  }
  points(com$cocoa$Z[1:10],v, type="l")
}

# Figure 3
setZDiscretization(10,0.7,com$cocoa)
lvls = c(5)
r = runif(200)
for (i in 1:200) {
  j=1
  while (r[i] > com$cocoa$Tz[lvls[i],j]) {
    r[i] = r[i]-com$cocoa$Tz[lvls[i],j]
    j = j+1
  }
  lvls[i+1] = j
}
plot(1:201, com$cocoa$Z[lvls], type="l")

# Set (Inverse) Demand function and real interest

setLinearPriceFunction = function(e) {
  e$P = function(x){return (e$a+e$b*x)}
  e$D = function(p){return ((p-e$a)/e$b)}
}

setRealInterest = function(r, e) {
  e$r = 0.05
}

for (e in com) {
  setLinearPriceFunction(e)
  setRealInterest(0.05, e)
}

# Transform to and from generic parameter vector theta using formula (43) (not the (43) that's actually (42))

toTheta = function(a,b,delta,e){
  e$a = a
  e$b = b
  e$delta = delta
  e$theta = c(a, log(-b), log(delta+0.05))
}

fromTheta = function(theta,e){
  e$theta = theta
  e$a = theta[1]
  e$b = -exp(theta[2])
  e$delta = -0.05+exp(theta[3])
}

initTheta = function (e) {
  toTheta(0.2, -0.15, 0.12, e)
}

for (e in com) {
  initTheta(e)
}

# 3.2 Discretization of the availability

setXDiscretization = function (m,e) {
  # minimum is 0 stored + minimum harvest, maximum is maximum harvest/delta with a safety factor allowing for delta to shrink
  safety_factor = 2
  e$X = (0:(m-1))/m*(e$Z[e$n]/e$delta*safety_factor-e$Z[1])+e$Z[1]
  e$m=m
}

setXDiscretization(20, com$cocoa)

# 3.2 Iterate to improve price function f
# f(x,z) is stored as a matrix with f[k,i] representing f(X[k],Z[i])

initF = function(e){
  e$f = matrix(nrow=e$m, ncol=e$n)
  for (k in 1:e$m){
    e$f[k,]=max(e$P(e$X[k]), 0)
  }
}


# Using formula (28)

iterateFOnce = function (e){
  newF = matrix(nrow=e$m, ncol=e$n)
  beta = (1-e$delta)/(1+e$r)
  splines=list()
  for (i in 1:e$n){
    splines=append(splines,splinefun(e$X, e$f[,i]))
  }
  futureValue = function(x,i){
    v = 0
    y = e$D(splines[[i]](x))
    for (j in 1:e$n){
      v = v + e$Tz[i,j]*splines[[j]](e$Z[j]+(1-e$delta)*(x-y))
    }
    return (beta*v)
  }
  for (k in 1:e$m){
    currentPrice=e$P(e$X[k])
    for (i in 1:e$n){
      newF[k,i]=max(futureValue(e$X[k],i),currentPrice)
    }
  }
  e$f=newF
}

iterateF = function(e, dCutoff = 1e-5, nCutoff = 100) {
  initF(e)
  lastF = e$f+Inf
  for (i in 1:nCutoff) {
    deviation = max(abs(e$f - lastF))
    if(deviation<dCutoff) {
      print(i)
      return (TRUE)
    }
    lastF = e$f
    iterateFOnce(e)
  }
  return (FALSE)
}


iterateF(com$cocoa)


plot(com$cocoa$X,com$cocoa$f[,1],type="l")
for(i in 2:com$cocoa$n){
  points(com$cocoa$X,com$cocoa$f[,i],type="l")
}



# 3.3 Generate the combined storage/harvest transition matrix

# calculates g(x,z) from formulas (36), (39)
expectedStorage = function(currentXLvl, currentZLvl, e) {
  return ((1-e$delta)*(e$X[currentXLvl]-e$D(e$f[currentXLvl,currentZLvl]))+e$rho*e$Z[currentZLvl])
}


# XZ Transition matrix will be ordered by Z on a large scale, X on a small scale
# X1Z1, X2Z1, ..., XmZ1, X1Z2, X2Z2, ..., xmZn

generateXZTransitionMatrix = function(e){
  theta = qnorm(0:e$n/e$n)
  e$Txz = matrix(nrow = e$m*e$n, ncol = e$m*e$n)
  xd = append(append(Inf,(e$X[2:e$m]-e$X[1:(e$m-1)])/2),Inf) # delta/2 between X, xd[i] is left of X[i], xd[i+1] right
  for (i in 1:e$m) {
    for (j in 1:e$n) {
      c_row = (i-1)*e$n + j
      for (k in 1:e$m) {
        for (l in 1:e$n) {
          c_column = (k-1)*e$n + l
          xdc = e$X[i] - expectedStorage(k,l,e)
          critXL = xdc - xd[i]
          critXR = xdc + xd[i+1]
          critZL = theta[j]-e$rho*e$Z[l]
          critZR = theta[j+1]-e$rho*e$Z[l]
          critL = max(critXL,critZL)
          critR = min(critXR,critZR)
          if (critL > critR){
            e$Txz[c_row,c_column]=0
          } else {
            e$Txz[c_row,c_column]=pnorm(critR)-pnorm(critL)
          }
        }
      }
    }
  }
}

for (e in com) {
  setZDiscretization(10,0.7,e)
  setXDiscretization(20,e)
  initTheta(e)
  iterateF(e)

  generateXZTransitionMatrix(e)
}


# Calculate invariant distribution

calculateInvariantDistribution = function(e) {
  e$id = eigen(e$Txz)$vectors[,1]
  if(Re(sum(e$id)) < 0) {
    e$id = - e$id
  }
}

for (e in com) {
  calculateInvariantDistribution(e)
}


# Visualize invariant distribution

s = com$cocoa

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

# 5.1 Estimate harvest level probabilities from price series

calculateGamma = function(e) {
  e$gamma = matrix(nrow=e$n, ncol=e$t)
  helper = matrix(nrow=e$n, ncol=e$t)    # contains non-normalized probabilities
  for (i in 1:e$n) {
    inversePriceFunction = splinefun(e$f[,i], e$X)
    probabilitySpline = splinefun(e$X, e$id[((i-1)*e$m+1):(i*e$m)])
    for (t in 1:e$t) {
      assumedX = inversePriceFunction(e$pDat[t])
      helper[i,t] = max(0, probabilitySpline(assumedX))
      # print(paste(i,t,assumedX,probabilitySpline(assumedX)))
    }
  }
  for (i in 1:e$n) {
    for (t in 1:e$t) {
      e$gamma[i,t] = helper[i,t]/sum(helper[,t])
    }
  }
  for (t in e$t){
    stopifnot(round(sum(e$gamma[,t]), digits=6) == 1)
  }
}


for (e in com){
  calculateGamma(e)
}



# 5.1 Calculate 1-period-ahead expectations and variances using formulas (45)/(46) and (48)/(49)


# Calculate Pseudo Log Likelihood Function


# Optimize Pseudo Log Likelihood Function

